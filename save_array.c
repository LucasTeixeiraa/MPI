#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "appctx.h"

void save_array(double *array, AppCtx *app)
{
    // save a parallel vtk XML image data format file
    int n    = app->global_n;
    int rank = app->rank;
    int size = app->size;
    int i, j;
    char filename[256];
    FILE *fp;

    if(rank != 0) return ;

    int x1 = 0;
    int x2 = n+1;
    int y1 = 0;
    int y2 = n+1;

    sprintf(filename, "array.%d.vti", rank);
    fp = fopen(filename, "w");
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", x1,x2,y1,y2);
    fprintf(fp, "<Piece Extent=\"%d %d %d %d 0 0\">\n", x1,x2,y1,y2);
    fprintf(fp, "<PointData>\n");
    // fprintf(fp, "<DataArray type=\"Float64\" Name=\"array\" NumberOfComponents=\"1\" format=\"binary\">\n");
    // fwrite(array, sizeof(double),app->ndof, fp);
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"array\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(int i=0;i<(n+2); i++) {
        for(int j=0;j <(n+2); j++)
            fprintf(fp, "%8.8e ", array[i*(n+2)+j]);
        fprintf(fp,"\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</ImageData>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);

} 
