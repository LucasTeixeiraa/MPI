
#include "appctx.h"

void save_array(double *array, AppCtx *app);


void ZeraVector(double *uold, double *unew, double *f, int ndof)
{
  for(int i=0;i<ndof;i++) {
    uold[i] = 0.0;
    unew[i] = 0.0;
    f[i]    = 0.0;
  }
}

// define uma macro para acessar os elementos do vetor
#define ind(i,j) (i)*(n+2)+(j)

int main(int argc, char **argv) 
{ 
  AppCtx app ;
  AppInit(&app,argc,argv);

  int n      = app.global_n;
  int energy = app.energy ;      // energy to be injected per iteration
  int niters = app.niters ;      // number of iterations

  int ndof   = app.ndof;         // number of degrees of freedom

  // Para a implementação paralela
  // Cada processo deve ter uma cópia local de uold e unew
  double *uold = (double*)malloc(ndof*sizeof(double)); // 
  double *unew = (double*)malloc(ndof*sizeof(double)); // 
  double *f    = (double*)malloc(ndof*sizeof(double)); // 
  double *tmp;

  // Inicializa os vetores
  ZeraVector(uold,unew,f,ndof);

  // Determina a posição dos pontos fonte
  int s1x = 1 + app.global_n / 4;
  int s1y = s1x;

  int s2x = 1 + 3 * app.global_n / 4;
  int s2y = s2x;

  // Injeta calor nos pontos fonte
  f[ind(s1x,s1y)] = app.energy;
  f[ind(s2x,s2y)] = -app.energy;

  // Cada processo deve determinar quais são os seus vizinhos
  // e alocar os buffers de envio e recebimento
  app.rank = 0;
  app.size = 1;
  app.neighbors[DOWN]  = -1;
  app.neighbors[UP]    = -1;
  app.neighbors[LEFT]  = -1;
  app.neighbors[RIGHT] = -1;
  app.coords[0] = 0;
  app.coords[1] = 0;
  
  double t=-AppGetTime();
 
  // Variaveris auxiliares
  double h = app.h;
  double h2= h*h;
  double error = 0.0;

  for(int iter=0; iter<niters; ++iter) 
  {

    // Cada processo deve trocar as bordas com seus vizinhos
    double local_error = 0.0;

    // Calculo do Stencil
    for(int i=1; i<= app.global_n; ++i) 
    {
      for(int j=1; j<= app.global_n; ++j) 
      {
        // 2nd order finite difference
        unew[ind(i,j)] = 0.25*(uold[ind(i-1,j)] + uold[ind(i+1,j)] + 
                               uold[ind(i,j-1)] + uold[ind(i,j+1)] + 
            
                                h2*f[ind(i,j)]);
        // calculo do erro local                     
        local_error += (unew[ind(i,j)]-uold[ind(i,j)])*(unew[ind(i,j)]-uold[ind(i,j)]);



      }
    }
    
    error = sqrt(local_error);

    printf("iter: %d   erro: %8.8e\n", iter,error);

    if(error < app.tol) break;

    tmp=unew; unew=uold; uold=tmp; // swap arrays
  }
  t+= AppGetTime();

  save_array(unew,&app);

  free(uold);
  free(unew);
  free(f);

  AppFinalize(&app);

  return 0;
}
