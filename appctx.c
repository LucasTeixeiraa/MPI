
#include "appctx.h"

/* Mostra a ajuda */
void show_help(char *name) {
    fprintf(stderr, "\
    [uso] %s <opcoes>\n\
        -h         mostra essa tela e sai.    \n\
        -n <int>   tamanho do grid.           \n\
        -i <int>   numero maximo de iterações \n\
        -t <float> tolerancia                 \n\
        -o <0|1|2> tipo de saida              \n\
            0: binario 1: ascii 2: VTI        \n\
        -e <int>   energia a ser injetado     .\n", name) ;
    exit(-1) ;
}

void AppInit(AppCtx* app, int argc, char*argv[])
{
   
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &app->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &app->size);
#else
    app->rank = 0;
    app->size = 1;
#endif
    AppGetParameters(app, argc, argv);
}


/* Inicializa o contexto da aplicação */
void AppGetParameters(AppCtx* app, int argc, char*argv[])
{
    // valores default
    app->global_n = 40;
    app->niters   = 1000;
    app->energy   = 10;
    app->L        = 1.0;
    app->output_type = OUTPUT_VTI;

    // processa as opções de linha de comando
    // modificando os valores defaults
    char opt;
    while( (opt = getopt(argc, argv, "hn:i:e:t:o:")) != -1 ) {
        switch ( opt ) {
            case 'h': /* help */
                show_help(argv[0]) ;
                break ;
            case 'n': /* opção -n */
                app->global_n = atoi(optarg) ;
                break ;
            case 'i': /* opção -i */
                app->niters = atoi(optarg) ;
                break ;
            case 'e': /* opção -e */
                app->energy = atoi(optarg) ;
                break ;
            case 't': /* opção -t */
                app->tol = atof(optarg) ;
                break ;
            case 'o': /* opção -o */
                app->output_type = atoi(optarg) ;
                if(app->output_type > 2 || app->output_type < 0) {
                    fprintf(stderr, "Tipo de saida invalida `%d'\n", app->output_type) ;
                    show_help(argv[0]) ;
                    AppFinalize(app);
                    exit(-1);
                }
                break ;
            default:
                fprintf(stderr, "Opcao invalida ou faltando argumento: `%c'\n", opt) ;
                AppFinalize(app);
                exit(-1) ;
        }
    }

    // calcula o tamanho local do grid
    app->local_n  = app->global_n/app->size;

    // calcula o número de graus de liberdade
    app->ndof     = (app->local_n+2)*(app->local_n+2);

    // calcula o tamanho do passo
    app->h        = app->L/(app->global_n+1);
}

void AppFinalize(AppCtx* app)
{
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

double AppGetTime()
{
#ifndef HAVE_MPI
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec*1e-6;
#else
    return MPI_Wtime();
#endif
}

