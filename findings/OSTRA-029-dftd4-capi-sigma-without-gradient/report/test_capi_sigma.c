#include <stdio.h>
#include <math.h>
#include "dftd4.h"
static double md(double*a,double*b,int n){double m=0;for(int i=0;i<n;i++){double d=fabs(a[i]-b[i]);if(d>m)m=d;}return m;}
int main(void){
  int num[4]={7,1,1,1};
  double xyz[12]={0,0,0.22, 0,1.77,-0.51, 1.53,-0.885,-0.51, -1.53,-0.885,-0.51};
  dftd4_error err=dftd4_new_error();
  dftd4_structure mol=dftd4_new_structure(err,4,num,xyz,NULL,NULL,NULL);
  dftd4_model disp=dftd4_new_d4_model(err,mol);
  dftd4_param par=dftd4_load_rational_damping(err,"pbe0",true);
  double e0,e2,g0[12],s0[9],s2[9];
  dftd4_get_dispersion(err,mol,disp,par,&e0,g0,s0);
  for(int i=0;i<9;i++) s2[i]=999.0;
  dftd4_get_dispersion(err,mol,disp,par,&e2,NULL,s2);
  printf("NH3: max|sigma|(both)=%.3e  max|sigma(sigma only)-sigma(both)|=%.3e  dE=%.1e\n",md(s0,(double[9]){0},9),md(s2,s0,9),e2-e0);
  return 0;}
