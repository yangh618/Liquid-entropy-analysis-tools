#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lattice.h"

int LoadLatticeFromXYZ(char * fname, lattice_t * lattice){
  char line[100];
  FILE *fid;

  if (!( fid = fopen(fname, "r"))) {
      printf("Can't find file : %s\n", fname);
      exit(0);
  }

  fscanf(fid, "%lf %lf %lf", &(lattice->ax), &(lattice->ay), &(lattice->az));
  fscanf(fid, "%lf %lf %lf", &(lattice->bx), &(lattice->by), &(lattice->bz));
  fscanf(fid, "%lf %lf %lf", &(lattice->cx), &(lattice->cy), &(lattice->cz));

  fscanf(fid, "%d", &(lattice->tot));
  fgets(line, 100, fid);
  
  int i;
  lattice->x = (double *) malloc(sizeof(double)*lattice->tot);
  lattice->y = (double *) malloc(sizeof(double)*lattice->tot);
  lattice->z = (double *) malloc(sizeof(double)*lattice->tot);
  for ( i=0;i<lattice->tot; i++) {
    fscanf(fid, "%lf %lf %lf", &(lattice->x[i]), &(lattice->y[i]), &(lattice->z[i]));
    fgets(line, 100, fid);
  }

  lattice->species = -1; 	/* no infos for specific elements */
  return 0;
}

int LoadLatticeFromPOSCAR(char *fname, lattice_t *lattice){
  char line[100];
  FILE *fid;

  if (!( fid = fopen(fname, "r"))) {
      printf("Can't find file : %s\n", fname);
      exit(0);
  }

  fgets(line, 100, fid);
  memcpy(lattice->comment, line, strlen(line));
  //fgets(line, 100, fid);
  double sfactor;
  fscanf(fid, "%lf", &sfactor);

  fscanf(fid, "%lf %lf %lf", &(lattice->ax), &(lattice->ay), &(lattice->az));
  fscanf(fid, "%lf %lf %lf", &(lattice->bx), &(lattice->by), &(lattice->bz));
  fscanf(fid, "%lf %lf %lf", &(lattice->cx), &(lattice->cy), &(lattice->cz));
  fgets(line, 100, fid);

  lattice->ax *=sfactor; lattice->ay *= sfactor; lattice->az *= sfactor;
  lattice->bx *=sfactor; lattice->by *= sfactor; lattice->bz *= sfactor;
  lattice->cx *=sfactor; lattice->cy *= sfactor; lattice->cz *= sfactor;
  
  lattice->species = 0;
  char c; int len;
  do{
    c = fgetc(fid);
    len=0;
    while ((( c >= 'A') && ( c<='Z')) || (( c >= 'a') && (c<='z'))){
      lattice->atom_sym[lattice->species][len] = c;
      len ++;
      c = fgetc(fid);
    }
    if ( len>0){
      lattice->atom_sym[lattice->species][len] = '\0';
      lattice->species ++ ;
    }
  }while ( c != '\n' );

  int i;
  lattice->tot = 0;
  for (i=0;i<lattice->species;i++) {
    fscanf(fid, "%d", &(lattice->atom_num[i]));
    lattice->tot += lattice->atom_num[i];
  }

  fgets(line, 100, fid);
  fgets(line, 100, fid);

  lattice->x = (double *) malloc(sizeof(double) *lattice->tot);
  lattice->y = (double *) malloc(sizeof(double) *lattice->tot);
  lattice->z = (double *) malloc(sizeof(double) *lattice->tot);
  for (i=0;i<lattice->tot; i++) fscanf(fid, "%lf %lf %lf", &(lattice->x[i]), &(lattice->y[i]), &(lattice->z[i]));
  
  return 0;
}

int LoadLatticeFromXDATCAR(char *fname, lattice_t *lattice){
  char line[100];
  FILE *fid;

  if (!( fid = fopen(fname, "r"))) {
      printf("Can't find file : %s\n", fname);
      exit(0);
  }

  fgets(line, 100, fid);
  memcpy(lattice->comment, line, strlen(line));
  //fgets(line, 100, fid);
  double sfactor;
  fscanf(fid, "%lf", &sfactor);

  fscanf(fid, "%lf %lf %lf", &(lattice->ax), &(lattice->ay), &(lattice->az));
  fscanf(fid, "%lf %lf %lf", &(lattice->bx), &(lattice->by), &(lattice->bz));
  fscanf(fid, "%lf %lf %lf", &(lattice->cx), &(lattice->cy), &(lattice->cz));
  fgets(line, 100, fid);

  lattice->ax *=sfactor; lattice->ay *= sfactor; lattice->az *= sfactor;
  lattice->bx *=sfactor; lattice->by *= sfactor; lattice->bz *= sfactor;
  lattice->cx *=sfactor; lattice->cy *= sfactor; lattice->cz *= sfactor;
  
  lattice->species = 0;
  char c; int len;
  do{
    c = fgetc(fid);
    len=0;
    while ((( c >= 'A') && ( c<='Z')) || (( c >= 'a') && (c<='z'))){
      lattice->atom_sym[lattice->species][len] = c;
      len ++;
      c = fgetc(fid);
    }
    if ( len>0){
      lattice->atom_sym[lattice->species][len] = '\0';
      lattice->species ++ ;
    }
  }while ( c != '\n' );

  int i;
  lattice->tot = 0;
  for (i=0;i<lattice->species;i++) {
    fscanf(fid, "%d", &(lattice->atom_num[i]));
    lattice->tot += lattice->atom_num[i];
  }

  fgets(line, 100, fid);
  fgets(line, 100, fid);

  lattice->x = (double *) malloc(sizeof(double) *lattice->tot);
  lattice->y = (double *) malloc(sizeof(double) *lattice->tot);
  lattice->z = (double *) malloc(sizeof(double) *lattice->tot);
  for (i=0;i<lattice->tot; i++) fscanf(fid, "%lf %lf %lf", &(lattice->x[i]), &(lattice->y[i]), &(lattice->z[i]));
  
  return 0;
}

int CopyLattice(lattice_t *dest, lattice_t *source){
  dest->ax = source->ax; dest->ay=source->ay; dest->az= source->az;
  dest->bx = source->bx; dest->by=source->by; dest->bz = source->bz;
  dest->cx = source->cx; dest->cy=source->cy; dest->cz = source->cz;

  int i;
  
  dest->tot = source->tot;
  if ( dest->x == NULL){
    dest->x = (double *) malloc(sizeof(double)*source->tot);
    dest->y = (double *) malloc(sizeof(double)*source->tot);
    dest->z = (double *) malloc(sizeof(double)*source->tot);
  }

  for (i=0;i<source->tot;i++) {
    dest->x[i] = source->x[i];
    dest->y[i] = source->y[i];
    dest->z[i] = source->z[i];
  }

  if (source->vx != NULL ){
    if ( dest->vx == NULL){
      dest->vx = (double *) malloc(sizeof(double)*source->tot);
      dest->vy = (double *) malloc(sizeof(double)*source->tot);
      dest->vz = (double *) malloc(sizeof(double)*source->tot);
    }
    
    for (i=0;i<source->tot;i++) {
      dest->vx[i] = source->vx[i];
      dest->vy[i] = source->vy[i];
      dest->vz[i] = source->vz[i];
    }
  }

  dest->species = source->species;
  for (i=0;i<source->species;i++) dest->atom_num[i] = source->atom_num[i];
  for (i=0;i<source->species;i++) dest->atom_radius[i] = source->atom_radius[i];
  for (i=0;i<source->species;i++) strcpy(dest->atom_sym[i], source->atom_sym[i]);

  return 0;
}

int PrintLattice(lattice_t *lattice){
  printf("%16.12lf %16.12lf %16.12lf\n", lattice->ax, lattice->ay, lattice->az);
  printf("%16.12lf %16.12lf %16.12lf\n", lattice->bx, lattice->by, lattice->bz);
  printf("%16.12lf %16.12lf %16.12lf\n", lattice->cx, lattice->cy, lattice->cz);

  int i;

  if ( lattice->species != -1){
    printf("%d\n", lattice->species);
    for (i=0;i<lattice->species;i++) printf("%s\t", lattice->atom_sym[i]);
    printf("\n");
    for (i=0;i<lattice->species;i++) printf("%d\t", lattice->atom_num[i]);
    printf("\n");
  }
  
  printf("%d\n", lattice->tot);
  for (i=0;i<lattice->tot;i++) printf("%16.12lf %16.12lf %16.12lf\n", lattice->x[i], lattice->y[i], lattice->z[i]);

  return 0;
}
int Print2Lattices(lattice_t *lattice, lattice_t *lattice2){
  printf("%16.12lf %16.12lf %16.12lf\t|\t %16.12lf %16.12lf %16.12lf\n", lattice->ax, lattice->ay, lattice->az, lattice2->ax, lattice2->ay, lattice2->az);
  printf("%16.12lf %16.12lf %16.12lf\t|\t %16.12lf %16.12lf %16.12lf\n", lattice->bx, lattice->by, lattice->bz, lattice2->bx, lattice2->by, lattice2->bz);
  printf("%16.12lf %16.12lf %16.12lf\t|\t %16.12lf %16.12lf %16.12lf\n", lattice->cx, lattice->cy, lattice->cz,lattice2->cx, lattice2->cy, lattice2->cz);

  int i;

  if ( lattice->species != -1){
    printf("%d\n", lattice->species);
    for (i=0;i<lattice->species;i++) printf("%s\t", lattice->atom_sym[i]);
    printf("\n");
    for (i=0;i<lattice->species;i++) printf("%d\t", lattice->atom_num[i]);
    printf("\n");
  }
  
  printf("%d\n", lattice->tot);
  for (i=0;i<lattice->tot;i++) printf("%d\t%16.12lf %16.12lf %16.12lf\t|\t %16.12lf %16.12lf %16.12lf\t%16.12lf\n", i+1, lattice->x[i], lattice->y[i], lattice->z[i], lattice2->x[i], lattice2->y[i], lattice2->z[i], sqrt(lattice2->x[i]*lattice2->x[i]+lattice2->y[i] *lattice2->y[i]+lattice2->z[i]*lattice2->z[i]));

  return 0;
}

int FreeLattice(lattice_t *lattice){
  free(lattice->x);
  free(lattice->y);
  free(lattice->z);
  free(lattice);
  return 0;
}

int LoadRWIGS(lattice_t * lattice){
  char *RWIGS="RWIGS";
  char *fname="POTCAR";
  char line[1000];
  char *head;
  double num;
  
  FILE *fid;
  
  if ( !(fid=fopen(fname,"r"))) printf("Can't find file %s\n", fname);

  int i,j,k;
  for (i=0;i<lattice->species;i++){
    while (!feof(fid)){
      fgets(line, 1000, fid);
      head = strstr(line, lattice->atom_sym[i]);
      if (head != NULL) break;
    }
    while (!feof(fid)){
      fgets(line, 1000, fid);
      head = strstr(line, RWIGS);
      if (head != NULL) break;
    }
    head = strstr(&(head[2]), RWIGS);
    k=0; num = 0;
    for (j=0;j<strlen(head);j++){
      if (head[j] == '.') k=1;
      if ((( head[j]>='0' ) && ( head[j]<='9' ))){
	num = num*10 + ((int)head[j]-(int)('0'));
	if (k>0) k++;
      }
    }
    if (k != 0) num=num/pow(10,k-1);
    lattice->atom_radius[i] = num;
  }
  return 0;
}

int LoadAtomMass(lattice_t * lattice){
  char *fMass="Mass";
  char *fname="POTCAR";
  char line[1000];
  char *head;
  double num;
  int ss;
  int i,j,k;
  
  FILE *fid;
  
  if ( (fid=fopen(fMass,"r"))) {
    fscanf(fid,"%d", &ss);
    for (i=0;i<lattice->species;i++) fscanf(fid,"%d", &ss);
    for (i=0;i<lattice->species;i++) fscanf(fid,"%lf", &(lattice->atom_mass[i]));
    fclose(fid);
    return 0;
  }
  if ( !(fid=fopen(fname,"r"))) {printf("Can't read Mass from POTCAR %s\n", fname); exit(0);}

  for (i=0;i<lattice->species;i++){
    while (!feof(fid)){
      fgets(line, 1000, fid);
      head = strstr(line, lattice->atom_sym[i]);
      if (head != NULL) break;
    }
    while (!feof(fid)){
      fgets(line, 1000, fid);
      head = strstr(line, "POMASS");
      if (head != NULL) break;
    }
    //head = strstr(&(head[2]), "POMASS");
    k=0; num = 0;
    for (j=0;j<strlen(head);j++){
      if (head[j] == '.') k=1;
      if ((( head[j]>='0' ) && ( head[j]<='9' ))){
	num = num*10 + ((int)head[j]-(int)('0'));
	if (k>0) k++;
      }
      if (head[j] == ';') break;
    }
    if (k != 0) num=num/pow(10,k-1);
    lattice->atom_mass[i] = num;
  }
  return 0;
}
