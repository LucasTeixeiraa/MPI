#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <stdbool.h>

#include "appctx.h"

/* Mostra a ajuda */
void show_help(char *name)
{
  fprintf(stderr, "\
            [uso] %s <opcoes>\n\
            -h         mostra essa tela e sai.    \n\
            -n <int>   tamanho do grid.           \n\
            -i <int>   numero maximo de iterações \n\
            -e <int>   energia a ser injetado     .\n",
          name);
  exit(-1);
}

void save_array(double *array, AppCtx *app);

void GetParameters(int argc, char *argv[], AppCtx *app)
{

  MPI_Comm_rank(MPI_COMM_WORLD, &app->rank);
  MPI_Comm_size(MPI_COMM_WORLD, &app->size);

  app->global_n = 300;
  app->niters = 1000;
  app->energy = 10;
  app->L = 1.0;
  char opt;
  char *optarg = NULL;
  while ((opt = getopt(argc, argv, "hn:i:e:")) != -1)
  {
    switch (opt)
    {
    case 'h': /* help */
      show_help(argv[0]);
      break;
    case 'n': /* opção -n */
      app->global_n = atoi(optarg);
      break;
    case 'i': /* opção -i */
      app->niters = atoi(optarg);
      break;
    case 'e': /* opção -e */
      app->energy = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Opcao invalida ou faltando argumento: `%c'\n", opt);
      exit(-1);
    }
  }

  app->local_n = app->global_n / app->size;

  app->ndof = (app->local_n + 2) * (app->global_n);
}

void ZeraVector(double *uold, double *unew, double *f, int ndof)
{
  for (int i = 0; i < ndof; i++)
  {
    uold[i] = 0.0;
    unew[i] = 0.0;
    f[i] = 0.0;
  }
}

// define uma macro para acessar os elementos do vetor

#define ind(i, j) (i) * (app.global_n) + (j)
#define ind_global(i, j) (i) * (app.global_n) + (j)

int main(int argc, char **argv)
{
  int tag1 = 1, tag2 = 1, tag3 = 1, tag4 = 1;
  int cont_com = 0;
  MPI_Init(&argc, &argv);

  AppCtx app;
  GetParameters(argc, argv, &app);

  FILE *file;
  char filename[20];
  sprintf(filename, "out%0d.txt", app.rank);
  file = fopen(filename, "w");

  MPI_Request request[4];
  MPI_Status stats[4];

  int n = app.global_n;    // número de pontos em cada dimensão
  int nl = app.local_n;    // número de pontos em cada dimensão que cada processo é responsável
  int energy = app.energy; // energy to be injected per iteration
  int niters = app.niters; // number of iterations

  int ndof = app.ndof; // number of degrees of freedom

  // Para a implementação paralela
  // Cada processo deve ter uma cópia local de uold e unew com tamanho (app->local_n + 2) * (app->global_n)
  double *uold = (double *)malloc(ndof * sizeof(double)); //
  double *unew = (double *)malloc(ndof * sizeof(double)); //
  double *f = (double *)malloc(ndof * sizeof(double));    //

  int buff_size = app.global_n; // tamanho do buffer usado para armazenar os valores

  double *recv_buff_down = (double *)malloc(buff_size * sizeof(double));
  double *recv_buff_up = (double *)malloc(buff_size * sizeof(double));

  double *send_buff_down = (double *)malloc(buff_size * sizeof(double));
  double *send_buff_up = (double *)malloc(buff_size * sizeof(double));

  double *tmp;

  double *usave = (double *)malloc((ndof * app.size) * sizeof(double)); // armazena a solução de todos os processos

  // Inicializa os vetores
  ZeraVector(uold, unew, f, ndof); // inicializa como zero
                                   //  Termo fonte
  int s1x = 1 + app.global_n / 4;
  int s1y = s1x;

  int s2x = 1 + 3 * app.global_n / 4;
  int s2y = s2x;
  
  // Cada processo deve determinar quais são os seus vizinhos
  // e alocar os buffers de envio e recebimento
  app.coords[0] = 0;        // todos tem a mesma coluna
  app.coords[1] = app.rank; // calcula a coordenda y  (linha)
  // Armazena os ranks dos vizinhos abaixo, acima, á esquerda e a direita. O -1 é usado para indicar a ausência de vizinho.
  app.neighbors[DOWN] = app.rank == app.size - 1 ? -1 : app.rank + 1; // Vizinho abaixo
  app.neighbors[UP] = app.rank == 0 ? -1 : app.rank - 1;              // Vizinho acima//Nesse cenário, não há necessidade de comunicação horizontal (esquerda/direita), pois todos os processos estão na mesma coluna global.

  // define quais indices da grade global o processo atual é responsável por atualizar
  int my_globX[2] = {0, n - 1};                                  // Limites das colunas englobadas
  int my_globY[2] = {app.rank * nl, app.rank * nl + nl - 1}; // Limites das linhas englobadas

  // verifica se os pontos fontes estão dentro do subdominio de responsabilidade do processo atual

  bool S1_here = s1y >= my_globY[0] && s1y <= my_globY[1];
  bool S2_here = s2y >= my_globY[0] && s2y <= my_globY[1];

  // ajustando os indices globais ao contexto local
  // subtraindo os limites inferiores dos indices globais
  if (S1_here)
  {
    printf("\nRank: %d - S1_here: %d\n", app.coords[1], S1_here);
    f[ind(s1x, s1y - my_globY[0])] = app.energy;
  }

  if (S2_here)
  {
    printf("\nRank: %d - S2_here: %d\n", app.coords[1], S2_here);
    f[ind(s2x, s2y - my_globY[0])] = -app.energy;
  }

  double heat = 0.0; // total heat in system

  double t = -MPI_Wtime();
  double sum_error = 0.0;

  double h = app.L / (double)(app.global_n + 1);
  double h2 = h * h;

  for (int iter = 0; iter < niters; ++iter)
  {
    // Definir os buffers de envio valores na borda do subdominio de cada processo que serão enviados ao processos vizinhos
    for (int i = 0; i < app.global_n; i++)
    {
      // Condições para processos que não são nem o primeiro nem o último
      if (app.rank != 0 && app.rank != app.size - 1)
      {
        send_buff_down[i] = uold[n * (nl) + i];
        send_buff_up[i] = uold[n + i];
      }
      // Condições para o primeiro processo
      else if (app.rank == 0)
      {
        send_buff_down[i] = uold[n * (nl) + i];
        // Para o processo no topo, pode definir o valor de acordo com a condição de contorno ou outra lógica
        send_buff_up[i] = -1; // ou outro valor conforme a necessidade
      }
      // Condições para o último processo
      else if (app.rank == app.size - 1)
      {
        // Para o processo na parte inferior, pode definir o valor de acordo com a condição de contorno ou outra lógica
        send_buff_down[i] = -1; // ou outro valor conforme a necessidade
        send_buff_up[i] = uold[n + i];
      }
    }


    cont_com = 0;
    MPI_Barrier(MPI_COMM_WORLD); // sincroniza os processos
                                 
    // Cada processo deve trocar as bordas com seus vizinhos
                                 
    MPI_Isend(send_buff_down, buff_size, MPI_DOUBLE, app.neighbors[DOWN], tag3, MPI_COMM_WORLD, &request[0]);

    MPI_Isend(send_buff_up, buff_size, MPI_DOUBLE, app.neighbors[UP], tag4, MPI_COMM_WORLD, &request[1]);

    MPI_Recv(recv_buff_down, buff_size, MPI_DOUBLE, app.neighbors[DOWN], tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(recv_buff_up, buff_size, MPI_DOUBLE, app.neighbors[UP], tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i < app.global_n; i++)
    {
      if (app.neighbors[DOWN] != -1)
        uold[n * (nl + 1) + i] = recv_buff_down[i];
      if (app.neighbors[UP] != -1)
        uold[i] = recv_buff_up[i];
    }
    
    double local_error = 0.0;

    // Calculo do Stencil
    for (int i = 1; i <= app.local_n; i++)
    {
      for (int j = 1; j < app.global_n - 1; j++) // Permanece global_n para as colunas
      {
        // 2nd order finite difference
        unew[ind(i, j)] = 0.25 * (uold[ind(i - 1, j)] + uold[ind(i + 1, j)] +
                                  uold[ind(i, j - 1)] + uold[ind(i, j + 1)] +
                                  h2 * f[ind(i, j)]);
        // calculo do erro local
        local_error += (unew[ind(i, j)] - uold[ind(i, j)]) * (unew[ind(i, j)] - uold[ind(i, j)]);
      }
    }

    MPI_Allreduce(&local_error, &sum_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    sum_error = sqrt(sum_error);

    if (app.rank == 0)
      printf("rank %d - iter: %d   sum_erro: %8.8e\n", app.rank, iter, sum_error);

    if (sum_error < 3.5e-6)
      break;

    tmp = unew;
    unew = uold;
    uold = tmp; // swap arrays
  }
  t += MPI_Wtime();

  fclose(file);

  if (app.rank == 0)
    printf("last heat: %f  erro: %8.8e  CPU time: %f\n", heat, sum_error, t);

  free(uold);
  free(unew);
  free(recv_buff_down);
  free(recv_buff_up);
  free(send_buff_down);
  free(send_buff_up);
  free(usave);
  MPI_Finalize();
  return 0;
}