/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  if (*kv) {
    for (int i = 0; i < *la; i++) {
      AB[4 * i] = 0;
      AB[4 * i + 1] = -1;
      AB[4 * i + 2] = 2;
      AB[4 * i + 3] = -1;
    }
  } else {
    for (int i = 0; i < *la; i++) {
      AB[3 * i] = -1;
      AB[3 * i + 1] = 2;
      AB[3 * i + 2] = -1;
    }
  }
  AB[*kv] = 0;
  AB[*la * (*lab - 1 + *kv) - 1] = 0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  for (size_t i = 0; i < *la; i++) {
    AB[1 + *kv] = 1;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int a;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (a=1;a<(*la)-1;a++){
    RHS[a]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int a;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (a=0;a<(*la);a++){
    EX_SOL[a] = (*BC0) + X[a]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int a;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (a=0;a<(*la);a++){
    x[a]=(a+1)*h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  return 0;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  double factor; 
  *info = 0;    
  if (*kl != 1 || *ku != 1) {
    printf("Erreur : La largeur des bandes doit être 1 pour une matrice tridiagonale.\n");
    *info = -1;
    return -1; 
  }
  for (int i = 0; i < *n - 1; i++) {
    if (AB[1 + i * (*lab)] == 0.0) { 
      *info = i + 1; 
      return -1;
    }
    factor = AB[0 + (i + 1) * (*lab)] / AB[1 + i * (*lab)];
    AB[0 + (i + 1) * (*lab)] = factor; 
    AB[1 + (i + 1) * (*lab)] -= factor * AB[2 + i * (*lab)];
  }
  for (int i = 0; i < *n; i++) {
    ipiv[i] = i + 1; 
  }
  return 0; 
}
