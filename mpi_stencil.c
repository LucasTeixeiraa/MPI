#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "appctx.h"

void save_array(double *array, AppCtx *app);

void ZeraVector(double *uold, double *unew, double *f, int ndof) {
  for (int i = 0; i < ndof; i++) {
    uold[i] = 0.0;
    unew[i] = 0.0;
    f[i] = 0.0;
  }
}

#define ind(i, j) ((i) * (n + 2) + (j))

int main(int argc, char **argv) {
  int tag1 = 1, tag2 = 2, tag3 = 3, tag4 = 4;
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  AppCtx app;
  AppInit(&app, argc, argv);

  int n = app.global_n;
  int energy = app.energy;
  int niters = app.niters;
  int ndof = app.ndof;
  int n_loc = app.local_n;

  int dims[2] = {0, 0};
  int periods[2] = {0, 0};
  MPI_Dims_create(size, 2, dims);
  MPI_Comm cart_comm;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);

  int coords[2];
  MPI_Cart_coords(cart_comm, rank, 2, coords);
  app.rank = rank;
  app.size = size;

  double *uold = (double *)malloc(ndof * sizeof(double));
  double *unew = (double *)malloc(ndof * sizeof(double));
  double *f_local = (double *)malloc(ndof * sizeof(double));
  double *tmp;
  ZeraVector(uold, unew, f_local, ndof);
  
  /*///////////////
    int buff_size = app.local_n + 2;

  double *recv_buff_left = (double*)malloc(buff_size*sizeof(double));
  double *recv_buff_right = (double*)malloc(buff_size*sizeof(double));
  double *recv_buff_down = (double*)malloc(buff_size*sizeof(double));
  double *recv_buff_up = (double*)malloc(buff_size*sizeof(double));

  double *send_buff_left = (double*)malloc(buff_size*sizeof(double));
  double *send_buff_right = (double*)malloc(buff_size*sizeof(double));
  double *send_buff_down = (double*)malloc(buff_size*sizeof(double));
  double *send_buff_up = (double*)malloc(buff_size*sizeof(double));
  ////////////////*/

  int s1x = 1 + n / 4;
  int s1y = 1 + n / 4;

  int s2x = 1 + 3 * n / 4;
  int s2y = 1 + 3 * n / 4;

  // Determine os pontos de início e fim de cada processo
  int x_start = coords[0] * n_loc + 1;
  int x_end = (coords[0] + 1) * n_loc;
  int y_start = coords[1] * n_loc + 1;
  int y_end = (coords[1] + 1) * n_loc;

  // Inserir os pontos de calor nos processos corretos
  if (s1x >= x_start && s1x <= x_end && s1y >= y_start && s1y <= y_end) {
    f_local[ind(s1x, s1y)] = app.energy;
  }
  if (s2x >= x_start && s2x <= x_end && s2y >= y_start && s2y <= y_end) {
    f_local[ind(s2x, s2y)] = -app.energy;
  }


  //double t = -AppGetTime();

  double h = app.h;
  //double h = app.L / ( double ) ( app.global_n + 1 );
  double sum_error = 0.0;
  double h2 = h * h;
  double error = 0.0;
  double temp=-MPI_Wtime();

  // Obter as coordenadas dos vizinhos
  int left, right, up, down;
  MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
  MPI_Cart_shift(cart_comm, 1, 1, &left, &right);

  for (int iter = 0; iter < niters; ++iter) {
    double local_error = 0.0;
    
    /*/////
        for(int i=1;i<=app.local_n;i++){
      send_buff_left[i] = uold[1+(app.local_n+2)*(i)];
      send_buff_right[i] = uold[(app.local_n)+(app.local_n+2)*(i)];  
      send_buff_down[i] = uold[(app.local_n + 2) + i];
      send_buff_up[i] = uold[(app.local_n+2)*app.local_n + i];  
    }
    //////////*/

    // Comunicar as bordas com os vizinhos
    if(app.neighbors[LEFT]!=-1) MPI_Sendrecv(&uold[ind(1, 0)], n + 2, MPI_DOUBLE, up, tag1,
                 &uold[ind(n + 1, 0)], n + 2, MPI_DOUBLE, down, tag1, cart_comm, MPI_STATUS_IGNORE);
    if(app.neighbors[RIGHT]!=-1) MPI_Sendrecv(&uold[ind(n, 0)], n + 2, MPI_DOUBLE, down, tag2,
                 &uold[ind(0, 0)], n + 2, MPI_DOUBLE, up, tag2, cart_comm, MPI_STATUS_IGNORE);
    if(app.neighbors[DOWN]!=-1) MPI_Sendrecv(&uold[ind(0, 1)], 1, MPI_DOUBLE, left, tag3,
                 &uold[ind(0, n + 1)], 1, MPI_DOUBLE, right, tag3, cart_comm, MPI_STATUS_IGNORE);
    if(app.neighbors[UP]!=-1) MPI_Sendrecv(&uold[ind(0, n)], 1, MPI_DOUBLE, right, tag4,
                 &uold[ind(0, 0)], 1, MPI_DOUBLE, left, tag4, cart_comm, MPI_STATUS_IGNORE);
    /*/////
      for(int i=1;i<=app.local_n;i++){
      uold[(app.local_n+2)*(i)] = recv_buff_left[i];
      uold[(app.local_n)+(app.local_n+2)*(i)+1] = recv_buff_right[i];  
      uold[i] = recv_buff_down[i];
      uold[(app.local_n+2)*(app.local_n + 1) + i] = recv_buff_up[i];  
    }
    /////////*/

    // Cálculo do Stencil
    for (int i = 1; i <= n; ++i) {
      for (int j = 1; j <= n; ++j) {
        unew[ind(i, j)] = 0.25 * (uold[ind(i - 1, j)] + uold[ind(i + 1, j)] +
                                  uold[ind(i, j - 1)] + uold[ind(i, j + 1)] +
                                  h2 * f_local[ind(i, j)]);
        local_error += (unew[ind(i, j)] - uold[ind(i, j)]) * (unew[ind(i, j)] - uold[ind(i, j)]);
      }
    }

   error = sqrt(local_error);

    
    MPI_Reduce(&error,&sum_error,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if (app.rank == 0) printf("rank %d - iter: %d   sum_erro: %8.8e\n", app.rank, iter, sum_error);

    //if (error < app.tol)
      //break;

    tmp = unew;
    unew = uold;
    uold = tmp;
  }

  // Sincronize todos os processos
  //MPI_Barrier(cart_comm);

  //t += AppGetTime();

  temp+=MPI_Wtime();
  save_array(unew, &app);
  if(app.rank == 0) printf("CPU time: %f\n",temp);
  free(uold);
  free(unew);
  free(f_local);

  AppFinalize(&app);

  MPI_Finalize();

  return 0;
}

