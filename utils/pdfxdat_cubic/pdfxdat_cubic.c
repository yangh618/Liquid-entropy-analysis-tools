#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef double vec_t[3];

int main(){
  int maxlen=300;
  char comment[maxlen];
  char *ptr, line[maxlen], endline[maxlen];

  fgets(comment, maxlen, stdin);
  double scale;
  scanf("%lf", &scale);

  double x,y,z,r,a, b, c;
  scanf("%lf %lf %lf", &x, &y, &z); a=x;
  scanf("%lf %lf %lf", &x, &y, &z); b=y;
  scanf("%lf %lf %lf", &x, &y, &z); c=z;
  fgets(endline, maxlen, stdin);

  int ntype=0, chem_numbers[100];
  char chem_symbols[100][3], *token=NULL;
  const char space[2] = " ";
  fgets(line, maxlen, stdin);
  line[strlen(line)-1]='\0';
  token = strtok(line, space);
  while ( token != NULL){
    if (( token != space )){
      strcpy(chem_symbols[ntype], token);
      ntype++;
    }
    token = strtok(NULL, space);
  }

  int i,j;
  /* for (i=0;i<ntype;i++) */
  /*   printf("%s\n", chem_symbols[i]); */

  int *chemhash, atomnum;
  atomnum=0;
  for (i=0;i<ntype;i++){
    scanf("%d", &(chem_numbers[i]));
    atomnum += chem_numbers[i];
  }
  chemhash = (int *) malloc(sizeof(int)*atomnum);
  atomnum=0;
  for (i=0;i<ntype;i++){
    for (j=atomnum; j< atomnum+chem_numbers[i]; j++)
      chemhash[j]=i;
    atomnum += chem_numbers[i];
  }
  vec_t *sites;

  int nsample=0;
  int nbin=2000;
  int **hist, t1, t2;

  hist = (int **) malloc(sizeof(int *)*(ntype+1)*ntype/2);
  for (i=0;i<(ntype+1)*ntype/2;i++)
    hist[i] = (int *) malloc(sizeof(int)* nbin);
  for (i=0;i<(ntype+1)*ntype/2;i++)
    for (j=0;j<nbin;j++)
      hist[i][j] = 0;
  
  sites = (vec_t *) malloc(sizeof(vec_t)*atomnum);
  while(fgets(line, maxlen, stdin) != NULL){
    if (*line=='\n')
      continue;
    printf("%s", line);
    if (line[0] == 'D'){
      for (i=0;i<atomnum;i++)
	scanf("%lf %lf %lf", &(sites[i][0]), &(sites[i][1]), &(sites[i][2]));
      //for (i=0;i<atomnum;i++)
      //printf("%lf %lf %lf\n", (sites[i][0]), (sites[i][1]), (sites[i][2]));
      for (i=0;i<atomnum;i++)
	for (j=0;j<i;j++){
	  x = sites[i][0]-sites[j][0];
	  y = sites[i][1]-sites[j][1];
	  z = sites[i][2]-sites[j][2];

	  if ( x > 0.5 ) x--; if (x<-0.5) x++;
	  if ( y > 0.5 ) y--; if (y<-0.5) y++;
	  if ( z > 0.5 ) z--; if (z<-0.5) z++;

	  r = sqrt(x*x*a*a +y*y*b*b +z*z*c*c);
	  //printf("%lf\n", floor(r*100));
	  t1 = chemhash[i] > chemhash[j] ? chemhash[i]:chemhash[j];
	  t2 = chemhash[i] < chemhash[j] ? chemhash[i]:chemhash[j];
	  
	  hist[t2*ntype+t1][(int)floor(r*100)]++;
	}
      nsample++;
    }
  }
  
  char fname[100]={'\0'};
  FILE *fid;
  double r2 ,r1, vol, dv, rho;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<=t1;t2++){
      sprintf(fname, "pdf.%s%s", chem_symbols[t2], chem_symbols[t1]);
      fid = fopen(fname, "w");
      vol = a*b*c;
      if (t1 == t2)
	rho = (double)chem_numbers[t1]*(chem_numbers[t1]-1)/2/vol/vol;
      else
	rho = (double)chem_numbers[t1]*(chem_numbers[t2])/vol;
      for (i=0;i<nbin;i++){
	r1 = (double)(i+0)/nbin*20;
	r2 = (double)(i+1)/nbin*20;
	dv = (r2*r2*r2-r1*r1*r1)/3*4*3.1415926;
	fprintf(fid, "%lf ", r2);
	fprintf(fid, "%lf\n", (double)hist[t2*ntype+t1][i]/nsample/vol/dv/rho);
      }
      fclose(fid);
    }
    
  free(chemhash);
}
