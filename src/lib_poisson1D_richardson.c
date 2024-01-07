/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
    int n= *la;
    int lda=n;
    double *AB = (double *)malloc(n * n * sizeof(double));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j) {
                AB[i * n + j] = 2.0;
            }else if(abs(i - j) == 1){
                AB[i * n + j] = -1.0;
            }else{
                AB[i * n + j] = 0.0;
        }
     }
    }
    int info;
    double work_query;
    int lwork = -1;
    LAPACK_dsyev("V", "U", &n, AB, &lda, eigval, &work_query, &lwork, &info);
    lwork = (int)work_query;
    double *work = (double *)malloc(lwork * sizeof(double));
    LAPACK_dsyev("V", "U", &n, AB, &lda, eigval, work, &lwork, &info);
    if(info != 0){
        fprintf(stderr, "Erreur dans la fonction LAPACK dsyev (info = %d)\n", info);
    }

    free(AB);
    free(work);
}


double eigmax_poisson1D(int *la) {
    int n = *la;
    int lda = n;
    double *AB = (double *)malloc(n * n * sizeof(double));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j) {
                AB[i * n + j] = 2.0;
            }else if(abs(i - j) == 1){
                AB[i * n + j] = -1.0;
            }else{
                AB[i * n + j] = 0.0;
        }
     }
    }

    int info;
    double work_query;
    int lwork = -1;
    LAPACK_dsyev("V", "U", &n, AB, &lda, NULL, &work_query, &lwork, &info);

    lwork = (int)work_query;
    double *work = (double *)malloc(lwork * sizeof(double));
    double *eigval = (double *)malloc(n * sizeof(double));

    LAPACK_dsyev("V", "U", &n, AB, &lda, eigval, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "Erreur dans la fonction LAPACK dsyev (info = %d)\n", info);
    }

    double eigmax = eigval[n - 1];

    free(AB);
    free(work);
    free(eigval);

    return eigmax;
}


double eigmin_poisson1D(int *la){
  return 0;
}

double richardson_alpha_opt(int *la){
  return 0;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

