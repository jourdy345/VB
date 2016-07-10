#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

double quadratur( double h, unsigned long int p, double q, double r, double eps, double del) //Funktion zur Berechnung der Startintegrale
{
  const int k_min=50; //Minimale Anzahl der Schritte bei der Quadratur //Minimale Anzahl der Schritte bei der Quadratur
  const int k_max=50000;
  double T, M, f, del2;
  int j=1;
  M=0.0;
  T=0.0;
  do
  {
    f=pow(j*h,p)*exp(q*j*h-r*j*j*h*h);
    T+=f;
    j++;
  }  while(j<k_min|| (f>del&& j<k_max));
  T=h*T;
  do
  {   M=pow(0.5*h,p)*exp(q*0.5*h-r*0.25*h*h);
          j=1;
    do
    {
      f=pow((j+0.5)*h,p)*exp(q*(j+0.5)*h-r*(j+0.5)*h*(j+0.5)*h);
      M+=f;
      j++;
    }while(j<k_min  || (f>del&& j<k_max));
    M=h*M;
    T=0.5*(T+M);
    h=0.5*h;
  }while(fabs(T-M)>eps);
return(T);
}
void R_delta( double* x, double* w, double* w2) {
 double q, r, h, I_0, I_1, quad, vorg, vorg2;
 double eps, del;
 unsigned long int p, k;
 vorg=0;
 p=x[0]; //Parameter
 q=x[1]; //Parameter
 r=x[2]; //Parameter
 h=1e-2; //Startschrittweite
 eps=1e-6; //Fehlerschranke |T-M|<eps
 del=1e-6; //Fehlerschranke "ab wann werden Folgeglieder der Trapezsumme vernachlässigt"
 //---------------------------Berechnung der Startintegrale------------------------------------
 //Berechnung von J(0,q,r)
I_0= quadratur( h,0,q,r,eps,del);
I_1= quadratur( h,1,q,r,eps,del);
 //---------------------------Iteration--------------------------------------------------------
 for (k=2;k<p+1;k++)
 { if(k==p-2) vorg2=I_1;
   if(k==p-1) vorg=I_1;
   quad=(q*I_1+(k-1)*I_0)/(2*r);
   I_0=I_1;
   I_1=quad;

 };
 
 w[0]=I_1/vorg;
  w2[0]=I_0/vorg;


}
// void foo() {}

R_CMethodDef cMethods[] = {
//  {"foo", (DL_FUNC) &foo, 1 }, {REALSXP} },
  {"R_delta", (DL_FUNC) &R_delta, 3 }, // {REALSXP,REALSXP,REALSXP} },
  {NULL, NULL, 0}
};

void R_init_VA(DllInfo *info)
{
/* Register the .C and .Call routines.
No .Fortran() or .External() routines,
so pass those arrays as NULL.
*/

R_registerRoutines(info,
cMethods,  NULL,
NULL, NULL);
}

void R_unload_VA(DllInfo *info)
{
  /* Release resources. */
}
