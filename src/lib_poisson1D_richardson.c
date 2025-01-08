/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
  double h = (1.0 / (double)*la + 1.0);
  for (int k = 0; k < *la; k++) {
    double sin_theta = sin((double)k * M_PI * h / 2.0);
    eigval[k-1] = 4.0 * sin_theta * sin_theta;
  }
}

double eigmax_poisson1D(int *la){
	double* eigval = (double*)malloc(*la * sizeof(double));
	eig_poisson1D(eigval, la);
	double eigmax = eigval[0];
	for (int i = 1; i < *la; i++) {
		if (eigval[i] > eigmax) {
			eigmax = eigval[i];
		}
	}
	return eigmax;
}

double eigmin_poisson1D(int *la){
	double* eigval = (double*)malloc(*la * sizeof(double));
	eig_poisson1D(eigval, la);
  double eigmin = eigval[0];
	for (int i = 1; i < *la; i++) {
		if (eigval[i] < eigmin) {
			eigmin = eigval[i];
		}
	}
	return eigmin;
}

double richardson_alpha_opt(int *la){
  return 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
	double* vec = malloc(sizeof(double) * *la);
	if (vec == NULL) {
		exit(EXIT_FAILURE);
	}

	memcpy(vec, RHS, *la * sizeof(double));
	cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, vec, 1);
	double norm_res = cblas_dnrm2(*la, vec, 1) / cblas_dnrm2(*la, RHS, 1);
	resvec[0] = norm_res;

	int iter = 1;
	for (iter = 1; iter < *maxit; ++iter) {
		cblas_daxpy(*la, *alpha_rich, vec, 1, X, 1);
		memcpy(vec, RHS, *la * sizeof(double));
		cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, vec, 1);

		norm_res = cblas_dnrm2(*la, vec, 1) / cblas_dnrm2(*la, RHS, 1);
		resvec[iter] = norm_res;

		if (norm_res < *tol) {
			break;
		}
	}
	*nbite = iter;
	free(vec);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
	X = calloc(*la, *la * sizeof(double));
	double* vec = malloc(sizeof(double) * *la);
	if (vec == NULL) {
		exit(EXIT_FAILURE);
	}

	memcpy(vec, RHS, *la * sizeof(double));
	cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, vec, 1);
	double norm_res = cblas_dnrm2(*la, vec, 1) / cblas_dnrm2(*la, RHS, 1);
	resvec[0] = norm_res;

	int iter = 1;
	for (iter = 1; iter < *maxit; ++iter) {
		cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, MB, *lab, vec, 1, 1.0, X, 1);
		memcpy(vec, RHS, *la * sizeof(double));
		cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, vec, 1);
		norm_res = cblas_dnrm2(*la, vec, 1) / cblas_dnrm2(*la, RHS, 1);
		resvec[iter] = norm_res;

		if (norm_res < *tol) {
			break;
		}
	}

	*nbite = iter;
	free(vec);
}

