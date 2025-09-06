#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cblas.h>

#include "lattice.h"

#define pi (double)3.1415926

double tol = 1e-2;
int fshort = 0;

typedef struct param_type{
  double *r;
  double *g;
  int len;
}param_t;

int GetFileSize(char *fname, int *len){
  FILE *fid;
  if (!(fid=fopen(fname, "r"))){
    printf("Cannot find file %s\n", fname);
    exit(0);
  }

  int length = 0;
  double x, y;
  while (!feof(fid)){
    fscanf(fid, "%lf %lf", &x, &y);
    //printf("%lf %lf\n", x, y);
    length++;
  }
  *len = length-1;
  //printf("%d\n", *len);

  return 0;
}

int Getgandr(char *fname, double *r, double *g){
  FILE *fid;
  if (!(fid=fopen(fname, "r"))){
    printf("Cannot find file %s\n", fname);
    exit(0);
  }

  int i =0 ;
  while (!feof(fid)){
    fscanf(fid, "%lf %lf", &(r[i]), &(g[i]));
    i++;
  }

  return 0;
}

int pdf2s2_myfunc(char *fname, int len, double *sfluct, double *sinfo, double *sum,
		  double rho, double alpha, double beta){
  double *r, *g, *logg;
  r = (double *) malloc(sizeof(double) * len);
  g = (double *) malloc(sizeof(double) * len);
  logg = (double *) malloc(sizeof(double) * len);

  FILE *fid;
  if (!(fid=fopen(fname, "r"))){
    printf("Cannot find file %s\n", fname);
    exit(0);
  }

  int i=0;
  while (!feof(fid)){
    fscanf(fid,"%lf %lf", &(r[i]), &(g[i]));
    i++;
  }
  fclose(fid);

  for ( i=0;i<len;i++ )
    if ( g[i] <= 0 ) logg[i] = 0;
    else logg[i] = log(g[i]);
  
  for ( i=0;i<len;i++ ){ sfluct[i] = 0; sinfo[i] = 0; sum[i] = 0;}

  for ( i=0;i<len-1;i++ ){
    sfluct[i+1] = sfluct[i] + 4*pi*r[i]*r[i]*((g[i]+g[i+1])/2-1)*(r[i+1]-r[i]);
    sinfo[i+1] = sinfo[i] + 4*pi*r[i]*r[i]*(-g[i]*logg[i]-g[i+1]*logg[i+1])/2*(r[i+1]-r[i]);
    sum[i+1] = sfluct[i+1] + sinfo[i+1];
  }

  for ( i=0;i<len-1;i++ ){
    sfluct[i+1] = 0.5*(alpha*beta + rho *alpha *beta * sfluct[i+1]);
    sinfo[i+1] = 0.5*rho *alpha *beta *sinfo[i+1];
    sum[i+1] = sfluct[i+1] + sinfo[i+1];
  }

  free(r); free(g); free(logg);
  return 0;
}

double fsfluct(double x, void *param_tmp){
  double val;
  
  param_t * param = (param_t *) param_tmp;

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *w = gsl_spline_alloc(gsl_interp_linear, param->len);
  gsl_spline_init(w, param->r, param->g, param->len);

  val = gsl_spline_eval(w, x, acc);
  gsl_spline_free(w);
  gsl_interp_accel_free(acc);

  return (val-1)*4*pi*x*x;
}

double fsinfo(double x, void *param_tmp){
  double val;
  
  param_t * param = (param_t *) param_tmp;

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *w = gsl_spline_alloc(gsl_interp_linear, param->len);
  gsl_spline_init(w, param->r, param->g, param->len);

  val = gsl_spline_eval(w, x, acc);
  gsl_spline_free(w);
  gsl_interp_accel_free(acc);

  if ( val <= 0 ) return 0;
  else return -val*log(val)*4*pi*x*x;
}

int gslintegration(double func(), int len, double *r, double *g, double r0, double *res, double *err){
  param_t * param = (param_t *) malloc(sizeof(param_t));
  param->r = r; param->g = g; param->len = len;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  gsl_function  F;
  F.function = func;
  F.params = param;

  gsl_integration_qag(&F, 0.01, r0, tol, tol, 1000, 6, w, res, err);
  
  gsl_integration_workspace_free(w);

  return 0;
}

int pdf2s2_gsl(int len, double *r, double *g, double r0, double rho, double alpha, double beta,
	       double *sfluct, double *sfluct_err, double *sinfo, double *sinfo_err,
	       double *sum, double *sum_err){

  gslintegration(fsfluct, len, r, g, r0, sfluct, sfluct_err);
  gslintegration(fsinfo, len, r, g, r0, sinfo, sinfo_err);
  //printf("%lf %lf %lf\n", r0, *sfluct, *sinfo);
  
  *sfluct = 0.5*(alpha*beta + rho *alpha *beta * *sfluct);
  *sfluct_err = 0.5* rho *alpha *beta * *sfluct_err;
  *sinfo = 0.5*rho *alpha *beta * *sinfo;
  *sinfo_err = 0.5*rho *alpha *beta * *sinfo_err;
  *sum = *sfluct + *sinfo;
  *sum_err = *sfluct_err + *sinfo_err;
  return 0;
}

int EntropyOfFILE(char *fname){
  int len;
  GetFileSize(fname, &len);
  
  return 0;
}

int EntropyOfOneFile_myfunc(char *fname){
  int k;

  double *sfluct, *sinfo, *sum, *r, *g;
  int len;
  double rho = 1;
  double alpha = 1, beta = 1;

  FILE *fid;
  
  GetFileSize(fname, &len);

  r = (double *) malloc(sizeof(double) * len);
  g = (double *) malloc(sizeof(double) * len);
  sfluct = (double *) malloc(sizeof(double) * len);
  sinfo = (double *) malloc(sizeof(double) * len);
  sum = (double *) malloc(sizeof(double) * len);
  
  Getgandr(fname, r, g);
  
  pdf2s2_myfunc(fname, len, sfluct, sinfo, sum, rho, alpha, beta);

  strcat(fname, ".s2");
  
  fid = fopen(fname, "w");
  for ( k=0;k<len;k++ ) fprintf(fid, "%lf %lf %lf %lf\n", r[k], sfluct[k], sinfo[k], sum[k]);
  fclose(fid);
  
  free(r); free(g);
  free(sfluct); free(sinfo); free(sum);
  return 0;
}

int EntropyOfAllFiles_myfunc(double temp, char *fxdat){
  int i,j,k;
  char fname[100];
  
  lattice_t *lattice=(lattice_t *) malloc(sizeof(lattice_t));
  LoadLatticeFromXDATCAR(fxdat, lattice);

  double *sfluct, *sinfo, *sum, *r, *g;
  double *sfluct_tot, *sinfo_tot, *sum_tot;
  int len;
  double vol = fabs(lattice->ax*lattice->by*lattice->cz-lattice->ax*lattice->bz*lattice->cy
    - lattice->ay*lattice->bx*lattice->cz+lattice->ay*lattice->bz*lattice->cx
		    + lattice->az*lattice->bx*lattice->cy-lattice->az*lattice->by*lattice->cx);
  double rho = ((double)lattice->tot) / vol;
  double alpha, beta, lambda;

  LoadAtomMass(lattice);
  double S1= 3.0/2;
  for (i=0;i<lattice->species;i++){
    alpha = ((double) lattice->atom_num[i]) / lattice->tot;
    //printf("%lf\n", lattice->atom_mass[i]);
    lambda = sqrt(1.239848*1.23984/(2*pi*lattice->atom_mass[i]*931.494e6*8.617e-5*temp))*1e4;
    S1 -= alpha*log(rho*alpha*lambda*lambda*lambda);
  }
  if ( fshort == 0 ){
    printf("Atomic Number: %d\n", lattice->tot);
    printf("Density [A^-3]: %lf\n", rho);
    printf("S1 (Sideal-1) [kB/atom]: %lf\n", S1);
  }
  else{
    printf("%d\n", lattice->tot);
    printf("%lf\n", rho);
    printf("%lf\n", S1);
  }
  FILE *fid;
  
  strcpy(fname, "pdf.");
  strcat(fname, lattice->atom_sym[0]);
  strcat(fname, lattice->atom_sym[0]);
  GetFileSize(fname, &len);
  sfluct_tot = (double *) malloc(sizeof(double) * len);
  sinfo_tot = (double *) malloc(sizeof(double) * len);
  sum_tot = (double *) malloc(sizeof(double) * len);
  for ( k=0;k<len;k++ ){sfluct_tot[k] = 0; sinfo_tot[k] = 0; sum_tot[k] = 0; }
  
  for ( i=0;i<lattice->species;i++ )
    for ( j=0;j<lattice->species;j++ ){
      strcpy(fname, "pdf.");
      strcat(fname, lattice->atom_sym[i]);
      strcat(fname, lattice->atom_sym[j]);
      //printf("%s\n", fname);
      GetFileSize(fname, &len);

      r = (double *) malloc(sizeof(double) * len);
      g = (double *) malloc(sizeof(double) * len);
      sfluct = (double *) malloc(sizeof(double) * len);
      sinfo = (double *) malloc(sizeof(double) * len);
      sum = (double *) malloc(sizeof(double) * len);

      Getgandr(fname, r, g);
      
      alpha = ((double)lattice->atom_num[i])/lattice->tot;
      beta = ((double)lattice->atom_num[j])/lattice->tot;
      pdf2s2_myfunc(fname, len, sfluct, sinfo, sum, rho, alpha, beta);

      strcpy(fname, "pdf.");
      strcat(fname, lattice->atom_sym[i]);
      strcat(fname, lattice->atom_sym[j]);
      strcat(fname, ".s2");

      fid = fopen(fname, "w");
      for ( k=0;k<len;k++ ) {
	fprintf(fid, "%lf %lf %lf %lf\n", r[k], sfluct[k], sinfo[k], sum[k]);
	sfluct_tot[k] += sfluct[k];
	sinfo_tot[k] += sinfo[k];
	sum_tot[k] += sum[k];
      }
      fclose(fid);

      free(r); free(g);
      free(sfluct); free(sinfo); free(sum);
    }

  strcpy(fname, "pdf.");
  strcat(fname, lattice->atom_sym[0]);
  strcat(fname, lattice->atom_sym[0]);
  r = (double *) malloc(sizeof(double) * len);
  g = (double *) malloc(sizeof(double) * len);
  Getgandr(fname, r, g);
  fid = fopen("tot.s2", "w");
  if ( lattice->species >= 1 ){
    for ( k=0;k<len;k++ ) 
      fprintf(fid, "%lf %lf %lf %lf\n", r[k], sfluct_tot[k], sinfo_tot[k], sum_tot[k]);
  }
  fclose(fid);
  free(r); free(g);
  free(sfluct_tot); free(sinfo_tot); free(sum_tot);
  
  return 0;
}

int EntropyOfOneFile_gsl(char *fname){
  int k;

  double *sfluct, *sinfo, *sum, *r, *g;
  int len;
  double rho = 1.0;
  double alpha = 1.0, beta = 1.0;

  FILE *fid;
  
  double *sfluct_err, *sinfo_err, *sum_err;
  int interv=100;
  GetFileSize(fname, &len);

  r = (double *) malloc(sizeof(double) * len);
  g = (double *) malloc(sizeof(double) * len);
  sfluct = (double *) malloc(sizeof(double) * interv);
  sinfo = (double *) malloc(sizeof(double) * interv);
  sum = (double *) malloc(sizeof(double) * interv);
  sfluct_err = (double *) malloc(sizeof(double) * interv);
  sinfo_err = (double *) malloc(sizeof(double) * interv);
  sum_err = (double *) malloc(sizeof(double) * interv);

  Getgandr(fname, r, g);
      
  for ( k=1;k<interv;k++){
    pdf2s2_gsl(len, r, g, r[len-1]/((double)interv)*k, rho, alpha, beta, &(sfluct[k]), &(sfluct_err[k]), &(sinfo[k]), &(sinfo_err[k]), &(sum[k]), &(sum_err[k]));
  }

  strcat(fname, ".s2");
  
  fid = fopen(fname, "w");
  for ( k=1;k<interv;k++ ) fprintf(fid, "%lf %lf %lf %lf %lf %lf %lf\n", r[len-1]/((double)interv)*k, sfluct[k], sinfo[k], sum[k], sfluct_err[k], sinfo_err[k], sum_err[k]);
  fclose(fid);
  
  free(r); free(g);
  free(sfluct); free(sinfo); free(sum);
  free(sfluct_err); free(sinfo_err); free(sum_err);

  return 0;
}


int EntropyOfAllFiles_gsl(char *fxdat){
  int i,j,k;
  char fname[100];
  
  lattice_t *lattice=(lattice_t *) malloc(sizeof(lattice_t));
  LoadLatticeFromXDATCAR(fxdat, lattice);

  double *sfluct, *sinfo, *sum, *r, *g;
  int len;
  double vol = fabs(lattice->ax*lattice->by*lattice->cz-lattice->ax*lattice->bz*lattice->cy
    - lattice->ay*lattice->bx*lattice->cz+lattice->ay*lattice->bz*lattice->cx
		    + lattice->az*lattice->bx*lattice->cy-lattice->az*lattice->by*lattice->cx);
  double rho = ((double)lattice->tot) / vol;
  double alpha, beta;

  FILE *fid;
  
  double *sfluct_err, *sinfo_err, *sum_err;
  int interv=100;

  double *sfluct_tot, *sinfo_tot, *sum_tot;
  double *sfluct_err_tot, *sinfo_err_tot, *sum_err_tot;
  
  strcpy(fname, "pdf.");
  strcat(fname, lattice->atom_sym[0]);
  strcat(fname, lattice->atom_sym[0]);
  GetFileSize(fname, &interv);
  sfluct_tot = (double *) malloc(sizeof(double) * interv);
  sinfo_tot = (double *) malloc(sizeof(double) * interv);
  sum_tot = (double *) malloc(sizeof(double) * interv);
  sfluct_err_tot = (double *) malloc(sizeof(double) * interv);
  sinfo_err_tot = (double *) malloc(sizeof(double) * interv);
  sum_err_tot = (double *) malloc(sizeof(double) * interv);
  for ( k=0;k<interv;k++ ){
    sfluct_tot[k] = 0; sinfo_tot[k] = 0; sum_tot[k] = 0;
    sfluct_err_tot[k] = 0; sinfo_err_tot[k] = 0; sum_err_tot[k] = 0;
  }

  for ( i=0;i<lattice->species;i++ )
    for ( j=0;j<lattice->species;j++ ){
      strcpy(fname, "pdf.");
      strcat(fname, lattice->atom_sym[i]);
      strcat(fname, lattice->atom_sym[j]);
      //printf("%s\n", fname);
      GetFileSize(fname, &len);

      r = (double *) malloc(sizeof(double) * len);
      g = (double *) malloc(sizeof(double) * len);
      sfluct = (double *) malloc(sizeof(double) * interv);
      sinfo = (double *) malloc(sizeof(double) * interv);
      sum = (double *) malloc(sizeof(double) * interv);
      sfluct_err = (double *) malloc(sizeof(double) * interv);
      sinfo_err = (double *) malloc(sizeof(double) * interv);
      sum_err = (double *) malloc(sizeof(double) * interv);

      Getgandr(fname, r, g);
      
      alpha = ((double)lattice->atom_num[i])/lattice->tot;
      beta = ((double)lattice->atom_num[j])/lattice->tot;

      for ( k=1;k<interv;k++){
	pdf2s2_gsl(len, r, g, r[len-1]/((double)interv)*k, rho, alpha, beta, &(sfluct[k]), &(sfluct_err[k]), &(sinfo[k]), &(sinfo_err[k]), &(sum[k]), &(sum_err[k]));
	sfluct_tot[k] += sfluct[k];
	sinfo_tot[k] += sinfo[k];
	sum_tot[k] += sum[k];
	sfluct_err_tot[k] += sfluct_err[k];
	sinfo_err_tot[k] += sinfo_err[k];
	sum_err_tot[k] += sum_err[k];
      }

      strcpy(fname, "pdf.");
      strcat(fname, lattice->atom_sym[i]);
      strcat(fname, lattice->atom_sym[j]);
      strcat(fname, ".s2");

      fid = fopen(fname, "w");
      for ( k=1;k<interv;k++ ) fprintf(fid, "%lf %lf %lf %lf %lf %lf %lf\n", r[len-1]/((double)interv)*k, sfluct[k], sinfo[k], sum[k], sfluct_err[k], sinfo_err[k], sum_err[k]);
      fclose(fid);

      free(r); free(g);
      free(sfluct); free(sinfo); free(sum);
      free(sfluct_err); free(sinfo_err); free(sum_err);
    }

  strcpy(fname, "pdf.");
  strcat(fname, lattice->atom_sym[0]);
  strcat(fname, lattice->atom_sym[0]);
  r = (double *) malloc(sizeof(double) * len);
  g = (double *) malloc(sizeof(double) * len);
  Getgandr(fname, r, g);
  fid = fopen("tot.s2", "w");
  if ( lattice->species >= 2 ){
    for ( k=1;k<interv;k++ )
      fprintf(fid, "%lf %lf %lf %lf %lf %lf %lf\n", r[k], sfluct_tot[k], sinfo_tot[k], sum_tot[k], sfluct_err_tot[k], sinfo_err_tot[k], sum_err_tot[k]);
  }
  fclose(fid);

  free(r); free(g);
  free(sfluct_tot); free(sinfo_tot); free(sum_tot);
  free(sfluct_err_tot); free(sinfo_err_tot); free(sum_err_tot);
  
  return 0;
}

int Usage(char * argv){
  printf("Calculate entropies for liquid states from two-body correlation functions.\n");

  printf("%s -x XDATCAR\n", argv);
  
  printf("\n+ Options:\n");
  printf("\t -h, --help: printf this page\n");
  printf("\t -x: XDATCAR\n");
  printf("\t -f: give a pdf file and compute its entropy assuming rho=1.0, alpha=1.0 and beta=1.0.\n");
  printf("\t -g: using gsl libs for integrations.\n");
  printf("\t -t: tolerance.\n ");
  
  printf("\n+ Input files:\n");
  printf("\t 1) XDATCAR. Atomic species and densities.\n");
  printf("\t 2) pdf files. pdf files should be named after two atomic symbols, eg. pdf.AlCo, pdf.AlCu and etc., and should go over all possible combinations of two arbitary elements.\n");
  printf("\t 3) Mass. This file include three lines. Line 1: the number of species. Line 2: atomic number for each species. Line 3: masses for each species\n");
  printf("\t 3) Trun. One number of temperature.\n");

  printf("\n+ Output files:\n");
  printf("\t stdout: S1 at temperature Trun\n");
  printf("\t pdf.XX.s2: fluctuation entropies, information entropies and total entropies (following by numerical errors from integration.\n");
  return 0;
}

int main(int argc, char **argv){
  int i;
  char fname[100];
  char fxdat[100];
  double temp;
  int file_only = 0;
  int gsl = 0;

  FILE *fid;
  if ( (fid=fopen("Trun", "r")) ){
    fscanf(fid, "%lf", &temp);
    fclose(fid);
  }
  
  /* command line options. */
  for ( i=1;i<argc;i++ ){
    if ((strcmp(argv[i], "-h")==0) || ( strcmp(argv[i], "--help")==0)) Usage(argv[0]);

    if ((strcmp(argv[i], "-f"))==0){
      i++;
      file_only = 1;
      strcpy(fname, argv[i]);
      continue;
    }

    if ((strcmp(argv[i], "-x"))==0){
      i++;
      strcpy(fxdat, argv[i]);
      continue;
    }

    if ((strcmp(argv[i], "-g")==0)){ gsl = 1 ;}

    if ((strcmp(argv[i], "-t")==0)){
      i++;
      tol = atof(argv[i]);
    }
    if ((strcmp(argv[i], "-T")==0)){ i++; temp = atof(argv[i]);}
    if ((strcmp(argv[i], "-s") == 0)) fshort=1;
  }
  
  if ((gsl == 1) && (file_only == 1)) EntropyOfOneFile_gsl(fname);
  if ((gsl == 1) && (file_only == 0)) EntropyOfAllFiles_gsl(fxdat);
  if ((gsl == 0) && (file_only == 1)) EntropyOfOneFile_myfunc(fname);
  if ((gsl == 0) && (file_only == 0)) EntropyOfAllFiles_myfunc(temp, fxdat);

  return 0;
}
