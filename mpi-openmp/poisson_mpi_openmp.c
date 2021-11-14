#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>

/*
 * Un paso del método de Jacobi para la ecuación de Poisson
 *
 *   Argumentos:
 *     - N,M: dimensiones de la malla
 *     - Entrada: x es el vector de la iteración anterior, b es la parte derecha del sistema
 *     - Salida: t es el nuevo vector
 *
 *   Se asume que x,b,t son de dimensión (N+2)*(M+2), se recorren solo los puntos interiores
 *   de la malla, y en los bordes están almacenadas las condiciones de frontera (por defecto 0).
 */
void jacobi_step(int N,int M,double *x,double *b,double *t, MPI_Comm *comm_cart)
{
  int ld = M+2;
  int rank;
  MPI_Comm_rank(*comm_cart, &rank);

  // Identificamos los vecinos de la malla

  enum DIRS {DOWN, UP, LEFT, RIGHT};
  int neighbours_ranks[4];

  MPI_Cart_shift( *comm_cart , 0 , 1 , &neighbours_ranks[LEFT] , &neighbours_ranks[RIGHT]);
  MPI_Cart_shift( *comm_cart , 1 , 1 , &neighbours_ranks[DOWN] , &neighbours_ranks[UP]);

  // Creamos el tipo para cuando mandemos columnas a la derecha e izquierda
  MPI_Datatype columna;
  MPI_Type_vector( N , 1 , ld , MPI_DOUBLE , &columna);
  MPI_Type_commit( &columna);

  // Envío de columnas
  if (rank%2 == 0){
    MPI_Send( &x[1*ld+M] , 1 , columna , neighbours_ranks[RIGHT] , 0 , *comm_cart);
    MPI_Recv( &x[1*ld+0] , 1 , columna , neighbours_ranks[LEFT] , 0 , *comm_cart , MPI_STATUS_IGNORE);
    MPI_Send( &x[1*ld+1] , 1 , columna , neighbours_ranks[LEFT] , 0 , *comm_cart);
    MPI_Recv( &x[1*ld+M+1] , 1 , columna , neighbours_ranks[RIGHT] , 0 , *comm_cart , MPI_STATUS_IGNORE);
  }
  else{
    MPI_Recv( &x[1*ld+0] , 1 , columna , neighbours_ranks[LEFT] , 0 , *comm_cart , MPI_STATUS_IGNORE);
    MPI_Send( &x[1*ld+M] , 1 , columna , neighbours_ranks[RIGHT] , 0 , *comm_cart);
    MPI_Recv( &x[1*ld+M+1] , 1 , columna , neighbours_ranks[RIGHT] , 0 , *comm_cart , MPI_STATUS_IGNORE);
    MPI_Send( &x[1*ld+1] , 1 , columna , neighbours_ranks[LEFT] , 0 , *comm_cart);
  }

  // Envío de filas
  if (rank%2 == 0){
    MPI_Send( &x[N*ld+1] , M , MPI_DOUBLE , neighbours_ranks[DOWN] , 0 , *comm_cart);
    MPI_Recv( &x[0*ld+1] , M , MPI_DOUBLE , neighbours_ranks[UP] , 0 , *comm_cart , MPI_STATUS_IGNORE);
    MPI_Send( &x[1*ld+1] , M , MPI_DOUBLE , neighbours_ranks[UP] , 0 , *comm_cart);
    MPI_Recv( &x[(N+1)*ld+1] , M , MPI_DOUBLE , neighbours_ranks[DOWN] , 0 , *comm_cart , MPI_STATUS_IGNORE);
  }
  else{
    MPI_Recv( &x[0*ld+1] , M , MPI_DOUBLE , neighbours_ranks[UP] , 0 , *comm_cart , MPI_STATUS_IGNORE);
    MPI_Send( &x[N*ld+1] , M , MPI_DOUBLE , neighbours_ranks[DOWN] , 0 , *comm_cart);
    MPI_Recv( &x[(N+1)*ld+1] , M , MPI_DOUBLE , neighbours_ranks[DOWN] , 0 , *comm_cart , MPI_STATUS_IGNORE);
    MPI_Send( &x[1*ld+1] , M , MPI_DOUBLE , neighbours_ranks[UP] , 0 , *comm_cart);
  }

  int i, j;
  #pragma omp parallel for private(j) collapse(2)
  for (i=1; i<=N; i++) {
    for (j=1; j<=M; j++) {
      t[i*ld+j] = (b[i*ld+j] + x[(i+1)*ld+j] + x[(i-1)*ld+j] + x[i*ld+(j+1)] + x[i*ld+(j-1)])/4.0;
    }
  }
}

/*
 * Método de Jacobi para la ecuación de Poisson
 *
 *   Suponemos definida una malla de (N+1)x(M+1) puntos, donde los puntos
 *   de la frontera tienen definida una condición de contorno.
 *
 *   Esta función resuelve el sistema Ax=b mediante el método iterativo
 *   estacionario de Jacobi. La matriz A no se almacena explícitamente y
 *   se aplica de forma implícita para cada punto de la malla. El vector
 *   x representa la solución de la ecuación de Poisson en cada uno de los
 *   puntos de la malla (incluyendo el contorno). El vector b es la parte
 *   derecha del sistema de ecuaciones, y contiene el término h^2*f.
 *
 *   Suponemos que las condiciones de contorno son igual a 0 en toda la
 *   frontera del dominio.
 */
void jacobi_poisson(int N,int M,double *x,double *b, MPI_Comm * comm_cart)
{
  int i, j, k, ld=M+2, conv, maxit=1000;
  double *t, local_s, total_s, tol=1e-6;

  t = (double*)calloc((N+2)*(M+2),sizeof(double));

  k = 0;
  conv = 0;

  int rank;
  MPI_Comm_rank(*comm_cart, &rank);

  while (!conv && k<maxit) {

    /* calcula siguiente vector */
    jacobi_step(N,M,x,b,t, comm_cart);

    /* criterio de parada: ||x_{k}-x_{k+1}||<tol */
    local_s = 0.0;

    #pragma omp parallel
    {
      #pragma omp for reduction(+:local_s) private(j) collapse(2)
      for (i=1; i<=N; i++) {
        for (j=1; j<=M; j++) {
          local_s += (x[i*ld+j]-t[i*ld+j])*(x[i*ld+j]-t[i*ld+j]);
        }
      }

      #pragma omp single nowait
      {
        MPI_Allreduce( &local_s , &total_s , 1 , MPI_DOUBLE , MPI_SUM , *comm_cart);
        conv = (sqrt(total_s)<tol);

        if (!rank){
          printf("Error en iteración %d: %g\n", k, sqrt(total_s));
        }

        k = k+1;
      }

      /* siguiente iteración */
      #pragma omp for private(j) collapse(2)
      for (i=1; i<=N; i++) {
        for (j=1; j<=M; j++) {
          x[i*ld+j] = t[i*ld+j];
        }
      }
    }
  }

  free(t);
}

int main(int argc, char **argv)
{
  int i, j, N=40, M=40, ld;
  double *x, *b, *sol, h=0.01, f=1.5;

  /* Extracción de argumentos */
  if (argc > 1) { /* El usuario ha indicado el valor de N */
    if ((N = atoi(argv[1])) < 0) N = 40;
  }
  if (argc > 2) { /* El usuario ha indicado el valor de M */
    if ((M = atoi(argv[2])) < 0) M = 1;
  }

  /* Se inician directivas MPI */
  MPI_Init( &argc , &argv);
  
  int size;
  MPI_Comm_size(MPI_COMM_WORLD , &size);

  int check_int = sqrt(size);
  float check_float = sqrt(size);

  int m,n;
  if (!(check_float - check_int)){
    m = M/check_int;
    n = N/check_int;
  }
  else{
    printf("La raíz del número de procesos ha de ser entera\n");
    exit(EXIT_FAILURE);
  }

  // Creación del comunicador cartesiano
  int dims[2] = {0,0};
  MPI_Dims_create( size , 2 , dims);
  int periods[2] = {0,0};
  int reorder = 1;

  MPI_Comm comm_cart;
  MPI_Cart_create( MPI_COMM_WORLD , 2 , dims , periods , reorder , &comm_cart);


  ld = m+2;  /* leading dimension */

  /* Reserva de memoria */
  x = (double*)calloc((n+2)*(m+2),sizeof(double));
  b = (double*)calloc((n+2)*(m+2),sizeof(double));

  /* Inicializar datos */
  //int max_threads = omp_get_max_threads();

  #pragma omp parallel for private(j) collapse(2)
  for (i=1; i<=n; i++) {
    for (j=1; j<=m; j++) {
      b[i*ld+j] = h*h*f;  /* suponemos que la función f es constante en todo el dominio */
    }
  }

  /* Resolución del sistema por el método de Jacobi */
  // Medimos tiempos de calculo de jacobi
  MPI_Barrier(comm_cart);
  double tic = MPI_Wtime();

  jacobi_poisson(n,m,x,b,&comm_cart);

  MPI_Barrier(comm_cart);
  double toc = MPI_Wtime();

  /* Recogida de la solución en máster */

  // Creamos tipo de dato bloque entero
  MPI_Datatype bloque, bloque_sol;
  MPI_Type_vector(n,m,ld,MPI_DOUBLE,&bloque);
  MPI_Type_commit(&bloque);
  MPI_Type_vector(n,m,M,MPI_DOUBLE,&bloque_sol);
  MPI_Type_commit(&bloque_sol);

  sol = (double*)calloc(N*M,sizeof(double));


  // identificamos rank, para definir el máster
  int rank;
  MPI_Comm_rank(comm_cart,&rank);

  int rank_cart;
  int coods_cart[2];
  /* Comunicación punto a punto */
  if (!rank){
    for (int i = 0; i<dims[0]; i++){
      for (int j = 0; j<dims[1]; j++){
        // Para la coordenada (0,0), es decirm rank = 0
        if (i == 0 && j == 0){
          for (int k=1; k<=n; k++){
            for (int p=1; p<=m; p++){
              sol[((dims[1]-1-j)*n*M + (k-1)*M) + p-1] = x[k*ld+p];
            }
          }
        }
        else{
          coods_cart[0] = i;
          coods_cart[1] = j;
          MPI_Cart_rank(comm_cart, coods_cart, &rank_cart);
          MPI_Recv(&sol[(dims[1]-1-j)*n*M + i*m],1,bloque_sol,rank_cart,0,comm_cart,MPI_STATUS_IGNORE);
        }

      }
    }
  }
  else{
    MPI_Send(&x[1*ld+1],1,bloque,0,0,comm_cart);
  }

  int num_threads;

  #pragma omp parallel
  {
    num_threads = omp_get_num_threads();
  }


  ld = M;

  /* Imprimir solución (solo para comprobación, eliminar en el caso de problemas grandes) */
  if (!rank){
    FILE *output = fopen("output.txt","w");

    fprintf(output, "Versión 'poisson.c' solo MPI\n");
    fprintf(output,"Tiempo de computo de la función 'jacobi_poisson': %f segundos\n", toc-tic);
    fprintf(output, "Tamaño: (N,M) = (%d, %d)\n", N, M);
    fprintf(output, "Número de threads usados: %d\n", num_threads);
    fprintf(output, "Número de procesos: %d\n",size);
    fclose(output);

    FILE *p = fopen("matrix_poisson_mpi_openmp.txt","w");
    if (!p)
      return 1;
    
    for (i=0; i<N; i++) {
      for (j=0; j<M; j++) {
        fprintf(p,"%g ", sol[i*ld+j]);
      }
      fprintf(p,"\n");
    }
    fclose(p);
  }
 

  MPI_Type_free(&bloque);
  MPI_Type_free(&bloque_sol);
  free(x);
  free(b);
  free(sol);

  MPI_Finalize();
  return 0;
}

