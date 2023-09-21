#ifndef APPCTX_H__
#define APPCTX_H__

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h>

 enum DIRECTIONS {DOWN=0, UP=1, LEFT=2, RIGHT=3};
 enum OUTPUT_TYPE {OUTPUT_BINARY=0, OUTPUT_ASCII=1, OUTPUT_VTI=2};

typedef struct 
{
    int global_n;      // global number of unknowns in each direction
    int local_n;       // local number of unknowns in each direction
    int ndof;          // number of degrees of freedom
    int niters;        // number of iterations
    int energy;        // energy to be injected per iteration
    int rank;          // rank do processo
    int size;           // numero de processos
    int neighbors[4];  // vizinhos
    int coords[2];     // coordenadas do processo
    double L;          // tamanho do dominio
    double h;
    double tol;        // tolerancia
    unsigned short     output_type; // tipo de saida        

} AppCtx; 

void AppInit(AppCtx *, int, char **);

void AppGetParameters(AppCtx *, int, char **);

void AppFinalize(AppCtx *);

double AppGetTime();

#endif /* APPCTX_H__ */
