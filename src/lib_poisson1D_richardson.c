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


double eigmin_poisson1D(int *la) {
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
    double eigmin = eigval[0];
    free(AB);
    free(work);
    free(eigval);

    return eigmin;
}


double richardson_alpha_opt(int *la) {
    int n = *la;
    int lda = n;
    
    double *AB = (double *)malloc(n * n * sizeof(double));
        for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j){
                AB[i * n + j]=2.0;
            }else if(abs(i - j)==1){
                AB[i * n + j]=-1.0;
            }else{
                AB[i * n + j]=0.0;
        }
     }
    }

    int info;
    double alpha_opt = 0.0;
    double *eigval = (double *)malloc(n * sizeof(double));
    double work_query;
    int lwork = -1;
    LAPACK_dsyev("V", "U", &n, AB, &lda, eigval, &work_query, &lwork, &info);
    lwork = (int)work_query;
    
    double *work = (double *)malloc(lwork * sizeof(double));
    LAPACK_dsyev("V", "U", &n, AB, &lda, eigval, work, &lwork, &info);
    alpha_opt = eigval[n - 1];
    
    free(AB);
    free(eigval);
    free(work);

    return alpha_opt;
}



void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    int i, j, it;
    double res, diff, normb;

    int N = *la + *ku + *kl;

    for (i = 0; i < N; i++) {
        X[i] = 0.0;
        resvec[i] = RHS[i];
        for (j = 0; j < N; j++){
            resvec[i] -= AB[i * (*lab) + j] * X[j];
    }
    }
    normb = 0.0;
    for (i = 0; i < N; i++){
        normb += RHS[i] * RHS[i];
        }
    normb = sqrt(normb);
    it = 0;
    do {
        for(i = 0; i < N; i++) {
            X[i] += *alpha_rich * resvec[i];
            resvec[i] = RHS[i];
            for(j = 0; j < N; j++){
                resvec[i] -= AB[i * (*lab) + j] * X[j];
                }
        }
        res = 0.0;
        for(i = 0; i < N; i++){
            res += resvec[i] * resvec[i];
            }
        res = sqrt(res);
        diff = 0.0;
        for(i = 0; i < N; i++){
            diff += fabs(*alpha_rich * resvec[i]);
            }
        it++;
        resvec[it] = res;

    }while(res/normb > *tol && it < *maxit);
    *nbite = it;
}


void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
    int N = *la + *ku + *kl;
    
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            MB[i*(*kv)+j]=0.0;
        }
    }
    
    for(int i=0;i<N;i++){
        MB[i*(*kv)+i]=AB[i*(*lab)+i];
        if(i<N-1){
            MB[i*(*kv)+i+1]=AB[i*(*lab)+i+1];
        }
        if(i>0){
            MB[i*(*kv)+i-1]=AB[i*(*lab)+i-1];
        }
    }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    int i,j;
    int N = *la + *ku + *kl;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            MB[i*(*kv)+j]=0.0;
            }
            }
            
    for(i=0;i<N;i++){
        for(j=0;j<=i;j++){
            MB[i*(*kv)+j]=AB[i*(*lab)+j];
        }
    }
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            MB[i*(*kv)+j]=AB[i*(*lab)+j];
        }
    }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    int i, j, it;
    double res, diff, normb;

    int N = *la + *ku + *kl;

    for (i = 0; i < N; i++){
        X[i] = 0.0;
	}
	
    for (i = 0; i < N; i++) {
        resvec[i] = RHS[i];
        for (j = 0; j < N; j++){
            resvec[i] -= AB[i * (*lab) + j] * X[j];
            }
    }
    
    normb = 0.0;
    for (i = 0; i < N; i++){
        normb += RHS[i] * RHS[i];
        }
    normb = sqrt(normb);
    it = 0;
    do {
        for (i = 0; i < N; i++){
            X[i] += resvec[i] / MB[i * (*la) + i];
            }
            
        for(i=0;i<N;i++){
            resvec[i]=RHS[i];
            for (j=0;j<N;j++){
                resvec[i] -= AB[i * (*lab) + j] * X[j];
                }
        }

        res = 0.0;
        for (i=0;i<N;i++){
            res +=resvec[i]*resvec[i];
            }
        res=sqrt(res);

        diff = 0.0;
        for(i=0;i<N;i++){
            diff += fabs(X[i]);
            }
        it++;
        resvec[it]=res;

    }while(res / normb > *tol && it < *maxit);
    *nbite = it;
}


