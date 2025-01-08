/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void plot_convergence_history(double* vec, int size, char* filename) {
  const char* prefix = "convergence_";
  const char* filename_full = filename + strlen(prefix);

  printf("Plotting the convergence history of %s method", filename_full);

  char data_filename[256];
  snprintf(data_filename, sizeof(data_filename), "%s_data.dat", filename);

  FILE* data_file = fopen(data_filename, "w");
  if (!data_file) {
    perror("Error opening data file");
    return;
  }

  for (int i = 0; i < size; i++) {
    fprintf(data_file, "%d %lf\n", i, vec[i]);
  }
  fclose(data_file);

  char gnuplot_script[10000];
  snprintf(gnuplot_script, sizeof(gnuplot_script),
    "set terminal png size 800,600;\n"
    "set output '%s.png';\n"
    "set title 'Convergence de la méthode de %s';\n"
    "set xlabel 'Itérations';\n"
    "set ylabel 'Erreur avant';\n"
    "plot '%s' using 1:2 with lines title 'Erreur';\n",
    filename,
    filename_full,
    data_filename);

  char script_filename[256];
  snprintf(script_filename, sizeof(script_filename), "%s_script.gnuplot", filename);

  FILE* script_file = fopen(script_filename, "w");
  if (!script_file) {
    perror("Error opening gnuplot script file");
    return;
  }
  fprintf(script_file, "%s", gnuplot_script);
  fclose(script_file);

  char command[512];
  snprintf(command, sizeof(command), "gnuplot %s", script_filename);
  int res = system(command);
  if (res) {
    perror("Error while running gnuplot");
  }

  remove(data_filename); 
  remove(script_filename); 
}