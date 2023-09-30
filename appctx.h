#ifndef APPCTX_H__
#define APPCTX_H__

enum DIRECTIONS {DOWN=0, UP=1, LEFT=2, RIGHT=3};

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
    //double tol;

} AppCtx; 

#endif /* APPCTX_H__ */