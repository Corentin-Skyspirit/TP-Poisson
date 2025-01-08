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
  return -2.0*cos((*la)*M_PI/(*la+1)) +2.0;
}

double eigmin_poisson1D(int *la){
  return -2.0*cos(M_PI/(*la+1)) +2.0;
}

double richardson_alpha_opt(int *la){
  return 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
	double* vec = malloc(*la * sizeof(double));
	if (vec == NULL) {
		exit(EXIT_FAILURE);
	}

	memcpy(vec, RHS, *la * sizeof(double));
	cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, vec, 1);
	double norm_res = cblas_dnrm2(*la, vec, 1) / cblas_dnrm2(*la, RHS, 1);
	resvec[0] = norm_res;

	int i = 1;
	for (i = 1; i < *maxit; ++i) {
		cblas_daxpy(*la, *alpha_rich, vec, 1, X, 1);
		memcpy(vec, RHS, *la * sizeof(double));
		cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, vec, 1);

		norm_res = cblas_dnrm2(*la, vec, 1) / cblas_dnrm2(*la, RHS, 1);
		resvec[i] = norm_res;

		if (norm_res < *tol) {
			break;
		}
	}
	*nbite = i;
	free(vec);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for (int i = 0; i < *la; i++) {
    for (int j = 0; j < *lab; j++) {
      int ind = (*lab)*i;
      MB[ind + j] = 0.0;
    }
    int ind = ((*lab)*i) + (*ku);
    MB[ind] = AB[ind];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for (int i=0; i<*la; i++){
    MB[*lab * i +1] = AB[i * (*lab) + 1];
    MB[*lab * i + 2] = AB[i*(*lab)+2];
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
	double* rk = malloc(sizeof(double) * (*la));
  double norm_B = cblas_dnrm2(*la, RHS, 1);
  double inv_norm = 1 / norm_B;
  double norm = 0.0;
  int* ipiv = malloc(sizeof(int) * (*la));
  int info = 0;
  int NHRS = 1;
  int kuu = (*ku)-1;

  dgbtrf_(la, la, kl, &kuu, MB, lab, ipiv, &info);
  for ((*nbite) = 0; (*nbite) < (*maxit); (*nbite)++) {

    for (int i = 0; i < (*la); i++) {
      rk[i] = RHS[i];
    }

    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, rk, 1);

    norm = cblas_dnrm2(*la, rk, 1) * inv_norm;
    resvec[(*nbite)] = norm;

    dgbtrs_("N", la, kl, &kuu, &NHRS, MB, lab, ipiv, rk, la, &info);

    cblas_daxpy(*la, 1, rk, 1, X, 1);

    if (norm <= (*tol))
      break;
  }

  free(rk);
  free(ipiv);
}

