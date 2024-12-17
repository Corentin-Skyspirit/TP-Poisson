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
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
