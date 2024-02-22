/** Code for external COMSOL material wrapper to ABAQUS umats:
    Small strain version
    Copyright (C) 2024  Sergio Lucarini */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

extern void umat_(double stress[6], double *statev, double ddsdde[6][6],
  double sse, double spd, double scd, double rpl, double ddsddt[6],
  double drplde[6], double drpldt, double stran[6], double dstran[6],
  double time[2], double *dtime, double temp, double dtemp, double predef[1],
  double dpred[1], char cmname, int *ndi, int *nshr, int *ntens, int *nstatv,
  double *props, int *nprops, double coords, double drot, double *pnewdt,
  double celent, double dfgrd0[3][3], double dfgrd1[3][3], int noel, int npt,
  int layer, int kspt, int kstep, int kinc);

EXPORT int eval(double e[6],         // Input: Strain tensor components in Voigt order (xx,yy,zz,yz,zx,xy)
                double s[6],         // Output: Cauchy stress components in Voigt order (xx,yy,zz,yz,zx,xy)
                double D[6][6],         // Output: Jacobian of stress with respect to strain, 6-by-6 matrix in row-major order
                int *nPar,         // Input: Number of material model parameters, scalar
                double *par,       // Input: Parameters
                int *nStates1,      // Input: Number of states, scalar        
                double *states1,     // States, nStates-vector
                char *errMsg // Error message argument, output, optional
                  ) { 		
  int ignoredi;
  double ignoredd;
  double ignoredd1[1];
  double ignoredd2[2];
  double ignoredd6[6];
  double ignoredd33[3][3];
  double time[2];
  double ddsdde[6][6];
  double stress[6];
  double dtime[1];
  int ndi = 3;
  int nshr = 3;
  int ntens = 6;
  int nstatev;
  int x,y,i,j,k,l,p,q,m;
  double pnewdt;
  double et[6];
  double de[6];
  double delta;

// initialise
  nstatev=nStates1[0]-13;
  pnewdt=1.0;
  time[0]=1.0*states1[nstatev];
  time[1]=1.0*states1[nstatev];
  dtime[0]=delta;
  
// strain in t and dstrain
  et[0]=1.0*states1[nstatev+1];
  et[1]=1.0*states1[nstatev+2];
  et[2]=1.0*states1[nstatev+3];
  et[3]=1.0*states1[nstatev+4];
  et[4]=1.0*states1[nstatev+5];
  et[5]=1.0*states1[nstatev+6];
  de[0]=e[0]-et[0];
  de[1]=e[1]-et[1];
  de[2]=e[2]-et[2];
  de[3]=(2*e[5]-et[3]);
  de[4]=(2*e[4]-et[4]);
  de[5]=(2*e[3]-et[5]);
  stress[0]=1.0*states1[nstatev+7];
  stress[1]=1.0*states1[nstatev+8];
  stress[2]=1.0*states1[nstatev+9];
  stress[3]=1.0*states1[nstatev+10];
  stress[4]=1.0*states1[nstatev+11];
  stress[5]=1.0*states1[nstatev+12];
  
// umat
    umat_(&stress[0],&states1[0],&ddsdde[0],ignoredd,ignoredd,ignoredd,ignoredd,
     ignoredd6,ignoredd6,ignoredd,&et[0],&de[0],&time[0],&dtime[0],
     ignoredd,ignoredd,ignoredd1,ignoredd1,ignoredd,
     &ndi,&nshr,&ntens,&nStates1[0],&par[0],&nPar[0],ignoredd,ignoredd,&pnewdt,ignoredd,
     ignoredd33,ignoredd33,ignoredi,ignoredi,ignoredi,ignoredi,
     ignoredi,ignoredi);
    if (pnewdt<1) {errMsg="no conv"; return 0;};
	
// save state
  states1[nstatev]=states1[10]+dtime[0];
  states1[nstatev+1]=1.*e[0];
  states1[nstatev+2]=1.*e[1];
  states1[nstatev+3]=1.*e[2];
  states1[nstatev+4]=2.*e[5];
  states1[nstatev+5]=2.*e[4];
  states1[nstatev+6]=2.*e[3];
  states1[nstatev+7]=1.0*stress[0];
  states1[nstatev+8]=1.0*stress[1];
  states1[nstatev+9]=1.0*stress[2];
  states1[nstatev+10]=1.0*stress[3];
  states1[nstatev+11]=1.0*stress[4];
  states1[nstatev+12]=1.0*stress[5];

// stress
  s[0]=1.0*stress[0];
  s[1]=1.0*stress[1];
  s[2]=1.0*stress[2];
  s[3]=1.0*stress[5];
  s[4]=1.0*stress[4];
  s[5]=1.0*stress[3];

// tangent
  D[0][0]=1.*ddsdde[0][0];D[0][1]=1.*ddsdde[0][1];D[0][2]=1.*ddsdde[0][2];D[0][3]=2.*ddsdde[0][5];D[0][4]=2.*ddsdde[0][4];D[0][5]=2.*ddsdde[0][3];
  D[1][0]=1.*ddsdde[1][0];D[1][1]=1.*ddsdde[1][1];D[1][2]=1.*ddsdde[1][2];D[1][3]=2.*ddsdde[1][5];D[1][4]=2.*ddsdde[1][4];D[1][5]=2.*ddsdde[1][3];
  D[2][0]=1.*ddsdde[2][0];D[2][1]=1.*ddsdde[2][1];D[2][2]=1.*ddsdde[2][2];D[2][3]=2.*ddsdde[2][5];D[2][4]=2.*ddsdde[2][4];D[2][5]=2.*ddsdde[2][3];
  D[3][0]=1.*ddsdde[5][0];D[3][1]=1.*ddsdde[5][1];D[3][2]=1.*ddsdde[5][2];D[3][3]=2.*ddsdde[5][5];D[3][4]=2.*ddsdde[5][4];D[3][5]=2.*ddsdde[5][3];
  D[4][0]=1.*ddsdde[4][0];D[4][1]=1.*ddsdde[4][1];D[4][2]=1.*ddsdde[4][2];D[4][3]=2.*ddsdde[4][5];D[4][4]=2.*ddsdde[4][4];D[4][5]=2.*ddsdde[4][3];
  D[5][0]=1.*ddsdde[3][0];D[5][1]=1.*ddsdde[3][1];D[5][2]=1.*ddsdde[3][2];D[5][3]=2.*ddsdde[3][5];D[5][4]=2.*ddsdde[3][4];D[5][5]=2.*ddsdde[3][3];
 
  return 0;  // Return value 0 if success, any other value trigger an exception
}
