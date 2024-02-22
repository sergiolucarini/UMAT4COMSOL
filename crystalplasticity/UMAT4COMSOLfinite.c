/** Code for external COMSOL material wrapper to ABAQUS umats:
    Finite strain version
    Copyright (C) 2023  Sergio Lucarini */

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
  double *props, int *nprops, double coords[3], double drot[3][3], double *pnewdt,
  double celent, double dfgrd0[3][3], double dfgrd1[3][3], int *noel, int npt,
  int layer, int kspt, int kstep, int kinc);

EXPORT int eval(double FlOld[3][3], // Deformation gradient at previous step, input
  double Fl[3][3], // Deformation gradient at current step, input
  double *tempOld, // Temperature at previous step, input
  double *temp, // Temperature at current step, input
  double *sysT, // Local material coordinate system, input
  double *delta, // Time step / Continutation parameter increment
  double Sl[6], // Second Piola-Kirchhoff stress, output
  double Jac[6][9], // Jacobian of stress with respect to deformation gradient,output
  int *nPar, // Number of material model parameters, input
  double *par, // Material model parameters, input
  int *nStates1, // Size of first state array, input, optional
  double *states1, // First state array, input/output, optional
  char *errMsg // Error message argument, output, optional
  ){
  int ignoredi;
  double ignoredd;
  double ignoredd1[1];
  double ignoredd2[2];
  double ignoredd6[6];
  double ignoredd33[3][3];
  double time[2];
  double ddsdde[6][6];
  double stress[6];
  double stran[6];
  double dstran[6];
  double drot[3][3];
  double dfgrd1[3][3];
  double dfgrd0[3][3];
  double coords[3];
  double dtime[1];
  int noel;
  double pnewdt;
  int ndi = 3;
  int nshr = 3;
  int ntens = 6;
  int x,y,i,j,k,l,p,q,m,n;
  int nstatev;
  double dsdf[3][3][3][3];
  double sum1[3][3][3][3];
  double sum11[3][3][3][3];
  double sum2[3][3][3][3];
  double sum3[3][3][3][3];
  double ddsdde4[3][3][3][3];
  double stre[3][3];
  double stra[3][3];
  double streS[3][3];
  double finv[3][3];
  double f2inv[3][3];
  double ftinv[3][3];
  double Indim[3][3];
  double det,detinv,dett,dettinv,detf2,detf2inv;
  double a2[3][3];
  double q2[3][3];
  double d2[3][3];
  double b2[3][3];
  double r2[3][3];
  double f2[3][3];
  double c2,s2,ssum,theta,determinant;

// initialise
  nstatev=nStates1[0]-22;
  pnewdt=1.0;
  noel=1*par[2];
  coords[0]=1.0*par[3];
  coords[1]=1.0*par[4];
  coords[2]=1.0*par[5];
  time[0]=1.0*states1[nstatev];
  time[1]=1.0*states1[nstatev];
  for (x = 0; x < 3; x++) {for (y = 0; y < 3; y++){dfgrd0[x][y]=1.0*states1[nstatev+1+3*x+y];};};
  dtime[0]=1.0*delta[0];
  if (time[0]==0){dfgrd0[0][0]=1.;dfgrd0[1][1]=1.;dfgrd0[2][2]=1.;};
  if (time[0]==0 && delta[0]==0){Fl[0][0]=1.;Fl[1][1]=1.;Fl[2][2]=1.;};
  dfgrd1[0][0]=1.*Fl[0][0];dfgrd1[0][1]=1.*Fl[1][0];dfgrd1[0][2]=1.*Fl[2][0];
  dfgrd1[1][0]=1.*Fl[0][1];dfgrd1[1][1]=1.*Fl[1][1];dfgrd1[1][2]=1.*Fl[2][1];
  dfgrd1[2][0]=1.*Fl[0][2];dfgrd1[2][1]=1.*Fl[1][2];dfgrd1[2][2]=1.*Fl[2][2];
  FlOld[0][0]=1.*dfgrd0[0][0];FlOld[0][1]=1.*dfgrd0[1][0];FlOld[0][2]=1.*dfgrd0[2][0];
  FlOld[1][0]=1.*dfgrd0[0][1];FlOld[1][1]=1.*dfgrd0[1][1];FlOld[1][2]=1.*dfgrd0[2][1];
  FlOld[2][0]=1.*dfgrd0[0][2];FlOld[2][1]=1.*dfgrd0[1][2];FlOld[2][2]=1.*dfgrd0[2][2];


// drot
  ftinv[0][0]=FlOld[1][1]*FlOld[2][2]-FlOld[2][1]*FlOld[1][2];
  ftinv[0][1]=FlOld[2][1]*FlOld[0][2]-FlOld[0][1]*FlOld[2][2];
  ftinv[0][2]=FlOld[0][1]*FlOld[1][2]-FlOld[1][1]*FlOld[0][2];
  ftinv[1][0]=FlOld[2][0]*FlOld[1][2]-FlOld[1][0]*FlOld[2][2];
  ftinv[1][1]=FlOld[0][0]*FlOld[2][2]-FlOld[2][0]*FlOld[0][2];
  ftinv[1][2]=FlOld[1][0]*FlOld[0][2]-FlOld[0][0]*FlOld[1][2];
  ftinv[2][0]=FlOld[1][0]*FlOld[2][1]-FlOld[2][0]*FlOld[1][1];
  ftinv[2][1]=FlOld[2][0]*FlOld[0][1]-FlOld[0][0]*FlOld[2][1];
  ftinv[2][2]=FlOld[0][0]*FlOld[1][1]-FlOld[1][0]*FlOld[0][1];
  dett=FlOld[0][0]*ftinv[0][0]+FlOld[0][1]*ftinv[1][0]+FlOld[0][2]*ftinv[2][0];
  for(m=0;m<3;m++){for(n=0;n<3;n++){ftinv[m][n]/=dett;};};
  for(i=0;i<3;i++){for (j=0;j<3;j++){f2[i][j]=0.0;
   for(k=0;k<3;k++){f2[i][j]+=Fl[i][k]*ftinv[k][j];};};};
  for(i=0;i<3;i++){for (j=0;j<3;j++){a2[i][j]=0.0;
   for(k=0;k<3;k++){a2[i][j]+=f2[k][i]*f2[k][j];};};};
  for(i=0;i<3;i++){for(j=0;j<3;j++){r2[i][j]=0.0;};r2[i][i]=1.0;};
  i=0;j=1;
  if(fabs(a2[0][2])<fabs(a2[0][1])){j=2;};
  ssum=1.0;
  while(ssum>0.0000000001){
     m=3-j;
     n=1-i;
     if(fabs(a2[0][m])>fabs(a2[n][2])){i=0;j=1*m;}
     else{i=1*n;j=2;};
     theta=0;
     if (a2[i][j]!= 0){theta=3.14159265358979323846/4.0;};
     if (a2[i][i]!=a2[j][j]){theta=0.5*atan(2.0*a2[i][j]/(a2[i][i]-a2[j][j]));};
     c2=cos(theta);s2=sin(theta);
     for(m=0;m<3;m++){for(n=0;n<3;n++){q2[m][n]=0.0;};q2[m][m]=1.0;};
     q2[i][i]=1.0*c2;q2[i][j]=1.0*s2;q2[j][i]=-1.0*s2;q2[j][j]=1.0*c2;
     for(m=0;m<3;m++){for(n=0;n<3;n++){d2[m][n]=0.0;
      for(k=0;k<3;k++){d2[m][n]+=q2[m][k]*a2[k][n];};};};
     for(m=0;m<3;m++){for(n=0;n<3;n++){a2[m][n]=0.0;
      for(k=0;k<3;k++){a2[m][n]+=d2[m][k]*q2[n][k];};};};
     for(m=0;m<3;m++){for(n=0;n<3;n++){b2[m][n]=0.0;
      for(k=0;k<3;k++){b2[m][n]+=q2[m][k]*r2[k][n];};};};
     for(m=0;m<3;m++){for(n=0;n<3;n++){r2[m][n]=1.0*b2[m][n];};};
     ssum=fabs(a2[0][1])+fabs(a2[0][2])+fabs(a2[1][2]);};
  a2[0][0]=sqrt(a2[0][0]);
  a2[1][1]=sqrt(a2[1][1]);
  a2[2][2]=sqrt(a2[2][2]);
  for(m=0;m<3;m++){for(n=0;n<3;n++){d2[m][n]=0.0;
   for(k=0;k<3;k++){d2[m][n]+=a2[m][k]*r2[k][n];};};};
  for(m=0;m<3;m++){for(n=0;n<3;n++){a2[m][n]=0.0;
   for(k=0;k<3;k++){a2[m][n]+=r2[k][m]*d2[k][n];};};};
  b2[0][0]=a2[1][1]*a2[2][2]-a2[2][1]*a2[1][2];
  b2[0][1]=a2[2][1]*a2[0][2]-a2[0][1]*a2[2][2];
  b2[0][2]=a2[0][1]*a2[1][2]-a2[1][1]*a2[0][2];
  b2[1][0]=a2[2][0]*a2[1][2]-a2[1][0]*a2[2][2];
  b2[1][1]=a2[0][0]*a2[2][2]-a2[2][0]*a2[0][2];
  b2[1][2]=a2[1][0]*a2[0][2]-a2[0][0]*a2[1][2];
  b2[2][0]=a2[1][0]*a2[2][1]-a2[2][0]*a2[1][1];
  b2[2][1]=a2[2][0]*a2[0][1]-a2[0][0]*a2[2][1];
  b2[2][2]=a2[0][0]*a2[1][1]-a2[1][0]*a2[0][1];
  determinant=a2[0][0]*b2[0][0]+a2[0][1]*b2[1][0]+a2[0][2]*b2[2][0];
  for(m=0;m<3;m++){for(n=0;n<3;n++){b2[m][n]/=determinant;};};
  for(m=0;m<3;m++){for(n=0;n<3;n++){r2[m][n]=0.0;
   for(k=0;k<3;k++){r2[m][n]+=f2[m][k]*b2[k][n];};};};
  drot[0][0]=1.0*r2[0][0];drot[0][1]=1.0*r2[1][0];drot[0][2]=1.0*r2[2][0];
  drot[1][0]=1.0*r2[0][1];drot[1][1]=1.0*r2[1][1];drot[1][2]=1.0*r2[2][1];
  drot[2][0]=1.0*r2[0][2];drot[2][1]=1.0*r2[1][2];drot[2][2]=1.0*r2[2][2];

// dstran
  for(x=0;x<3;x++){for(y=0;y<3;y++){f2[x][y]=(Fl[x][y]+FlOld[x][y])/2.0;};};
  f2inv[0][0]=f2[1][1]*f2[2][2]-f2[2][1]*f2[1][2];
  f2inv[0][1]=f2[2][1]*f2[0][2]-f2[0][1]*f2[2][2];
  f2inv[0][2]=f2[0][1]*f2[1][2]-f2[1][1]*f2[0][2];
  f2inv[1][0]=f2[2][0]*f2[1][2]-f2[1][0]*f2[2][2];
  f2inv[1][1]=f2[0][0]*f2[2][2]-f2[2][0]*f2[0][2];
  f2inv[1][2]=f2[1][0]*f2[0][2]-f2[0][0]*f2[1][2];
  f2inv[2][0]=f2[1][0]*f2[2][1]-f2[2][0]*f2[1][1];
  f2inv[2][1]=f2[2][0]*f2[0][1]-f2[0][0]*f2[2][1];
  f2inv[2][2]=f2[0][0]*f2[1][1]-f2[1][0]*f2[0][1];
  detf2=f2[0][0]*f2inv[0][0]+f2[0][1]*f2inv[1][0]+f2[0][2]*f2inv[2][0];
  for(m=0;m<3;m++){for(n=0;n<3;n++){f2inv[m][n]/=detf2;};};
  for(x = 0; x < 6; x++){dstran[x]=0.;};
  for(k = 0; k < 3; k++){
     dstran[0]+=(Fl[0][k]-FlOld[0][k])*f2inv[k][0];
     dstran[1]+=(Fl[1][k]-FlOld[1][k])*f2inv[k][1];
     dstran[2]+=(Fl[2][k]-FlOld[2][k])*f2inv[k][2];
     dstran[3]+=(Fl[0][k]-FlOld[0][k])*f2inv[k][1]+(Fl[1][k]-FlOld[1][k])*f2inv[k][0];
     dstran[4]+=(Fl[0][k]-FlOld[0][k])*f2inv[k][2]+(Fl[2][k]-FlOld[2][k])*f2inv[k][0];
     dstran[5]+=(Fl[1][k]-FlOld[1][k])*f2inv[k][2]+(Fl[2][k]-FlOld[2][k])*f2inv[k][1];};

// stran and stress in t
  stre[0][0]=1.0*states1[nstatev+10];
  stre[1][1]=1.0*states1[nstatev+11];
  stre[2][2]=1.0*states1[nstatev+12];
  stre[0][1]=1.0*states1[nstatev+13];
  stre[1][0]=1.0*states1[nstatev+13];
  stre[0][2]=1.0*states1[nstatev+14];
  stre[2][0]=1.0*states1[nstatev+14];
  stre[1][2]=1.0*states1[nstatev+15];
  stre[2][1]=1.0*states1[nstatev+15];
  stra[0][0]=1.0*states1[nstatev+16];
  stra[1][1]=1.0*states1[nstatev+17];
  stra[2][2]=1.0*states1[nstatev+18];
  stra[0][1]=0.5*states1[nstatev+19];
  stra[1][0]=0.5*states1[nstatev+19];
  stra[0][2]=0.5*states1[nstatev+20];
  stra[2][0]=0.5*states1[nstatev+20];
  stra[1][2]=0.5*states1[nstatev+21];
  stra[2][1]=0.5*states1[nstatev+21];
  for(x = 0; x < 6; x++){stress[x]=0.0;};
  for(x = 0; x < 6; x++){stran[x]=0.0;};
  for(k = 0; k < 3; k++){
     for (l = 0; l < 3; l++){
         stress[0]+=r2[0][k]*stre[k][l]*r2[0][l];
         stress[1]+=r2[1][k]*stre[k][l]*r2[1][l];
         stress[2]+=r2[2][k]*stre[k][l]*r2[2][l];
         stress[3]+=0.5*(r2[0][k]*stre[k][l]*r2[1][l]+r2[1][k]*stre[k][l]*r2[0][l]);
         stress[4]+=0.5*(r2[0][k]*stre[k][l]*r2[2][l]+r2[2][k]*stre[k][l]*r2[0][l]);
         stress[5]+=0.5*(r2[1][k]*stre[k][l]*r2[2][l]+r2[2][k]*stre[k][l]*r2[1][l]);
         stran[0]+=r2[0][k]*stra[k][l]*r2[0][l];
         stran[1]+=r2[1][k]*stra[k][l]*r2[1][l];
         stran[2]+=r2[2][k]*stra[k][l]*r2[2][l];
         stran[3]+=r2[0][k]*stra[k][l]*r2[1][l]+r2[1][k]*stra[k][l]*r2[0][l];
         stran[4]+=r2[0][k]*stra[k][l]*r2[2][l]+r2[2][k]*stra[k][l]*r2[0][l];
         stran[5]+=r2[1][k]*stra[k][l]*r2[2][l]+r2[2][k]*stra[k][l]*r2[1][l];};};

// umat
     umat_(&stress[0],&states1[0],&ddsdde[0],ignoredd,ignoredd,ignoredd,ignoredd,
      ignoredd6,ignoredd6,ignoredd,&stran[0],&dstran[0],&time[0],&dtime[0],
      ignoredd,ignoredd,ignoredd1,ignoredd1,ignoredd,
      &ndi,&nshr,&ntens,&nStates1[0],&par[0],&nPar[0],&coords[0],&drot[0],&pnewdt,ignoredd,
      &dfgrd0[0],&dfgrd1[0],&noel,ignoredi,ignoredi,ignoredi,
      ignoredi,ignoredi);
      if (pnewdt<1) {errMsg="no conv"; return 0;}; 

// save state
  states1[nstatev]=states1[nstatev]+dtime[0];
  for (x = 0; x < 3; x++) {for (y = 0; y < 3; y++){states1[nstatev+1+3*x+y]=1.0*dfgrd1[x][y];}};
  for (x = 0; x < 6; x++) {states1[nstatev+10+x]=1.0*stress[x];};
  for (x = 0; x < 6; x++) {states1[nstatev+16+x]=stran[x]+dstran[x];};

// stress
  stre[0][0]=1.0*stress[0];
  stre[1][1]=1.0*stress[1];
  stre[2][2]=1.0*stress[2];
  stre[0][1]=1.0*stress[3];
  stre[1][0]=1.0*stress[3];
  stre[0][2]=1.0*stress[4];
  stre[2][0]=1.0*stress[4];
  stre[1][2]=1.0*stress[5];
  stre[2][1]=1.0*stress[5];
  finv[0][0]=detinv*(Fl[1][1]*Fl[2][2]-Fl[1][2]*Fl[2][1]);
  finv[1][1]=detinv*(Fl[0][0]*Fl[2][2]-Fl[0][2]*Fl[2][0]);
  finv[2][2]=detinv*(Fl[0][0]*Fl[1][1]-Fl[0][1]*Fl[1][0]);
  finv[0][1]=-detinv*(Fl[0][1]*Fl[2][2]-Fl[2][1]*Fl[0][2]);
  finv[1][0]=-detinv*(Fl[1][0]*Fl[2][2]-Fl[1][2]*Fl[2][0]);
  finv[0][2]=detinv*(Fl[0][1]*Fl[1][2]-Fl[0][2]*Fl[1][1]);
  finv[2][0]=detinv*(Fl[1][0]*Fl[2][1]-Fl[1][1]*Fl[2][0]);
  finv[1][2]=-detinv*(Fl[0][0]*Fl[1][2]-Fl[0][2]*Fl[1][0]);
  finv[2][1]=-detinv*(Fl[0][0]*Fl[2][1]-Fl[0][1]*Fl[2][0]);
  finv[0][0]=Fl[1][1]*Fl[2][2]-Fl[2][1]*Fl[1][2];
  finv[0][1]=Fl[2][1]*Fl[0][2]-Fl[0][1]*Fl[2][2];
  finv[0][2]=Fl[0][1]*Fl[1][2]-Fl[1][1]*Fl[0][2];
  finv[1][0]=Fl[2][0]*Fl[1][2]-Fl[1][0]*Fl[2][2];
  finv[1][1]=Fl[0][0]*Fl[2][2]-Fl[2][0]*Fl[0][2];
  finv[1][2]=Fl[1][0]*Fl[0][2]-Fl[0][0]*Fl[1][2];
  finv[2][0]=Fl[1][0]*Fl[2][1]-Fl[2][0]*Fl[1][1];
  finv[2][1]=Fl[2][0]*Fl[0][1]-Fl[0][0]*Fl[2][1];
  finv[2][2]=Fl[0][0]*Fl[1][1]-Fl[1][0]*Fl[0][1];
  det=Fl[0][0]*finv[0][0]+Fl[0][1]*finv[1][0]+Fl[0][2]*finv[2][0];
  for(m=0;m<3;m++){for(n=0;n<3;n++){finv[m][n]/=det;};};
  for (i = 0; i < 3; i++) {
     for (j = 0; j < 3; j++) {
         streS[i][j]=0;
         for (q = 0; q < 3; q++){
             for (p = 0; p < 3; p++){
                 streS[i][j]+=finv[i][p]*stre[p][q]*finv[j][q];};};
         streS[i][j]=det*streS[i][j];};};
  Sl[0]=1.0*streS[0][0];
  Sl[1]=1.0*streS[1][1];
  Sl[2]=1.0*streS[2][2];
  Sl[3]=(streS[1][2]+streS[2][1])/2.0;
  Sl[4]=(streS[0][2]+streS[2][0])/2.0;
  Sl[5]=(streS[0][1]+streS[1][0])/2.0;

// tangent
  Indim[0][0]=1.0;
  Indim[0][1]=0.0;
  Indim[0][2]=0.0;
  Indim[1][0]=0.0;
  Indim[1][1]=1.0;
  Indim[1][2]=0.0;
  Indim[2][0]=0.0;
  Indim[2][1]=0.0;
  Indim[2][2]=1.0;
  ddsdde4[0][0][0][0]=1.0*ddsdde[0][0];
  ddsdde4[0][0][1][1]=1.0*ddsdde[1][0];
  ddsdde4[1][1][0][0]=1.0*ddsdde[0][1];
  ddsdde4[1][1][1][1]=1.0*ddsdde[1][1];
  ddsdde4[0][0][2][2]=1.0*ddsdde[2][0];
  ddsdde4[2][2][0][0]=1.0*ddsdde[0][2];
  ddsdde4[1][1][2][2]=1.0*ddsdde[2][1];
  ddsdde4[2][2][1][1]=1.0*ddsdde[1][2];
  ddsdde4[2][2][2][2]=1.0*ddsdde[2][2];
  ddsdde4[0][0][0][1]=1.0*ddsdde[3][0];
  ddsdde4[0][0][1][0]=1.0*ddsdde[3][0];
  ddsdde4[0][0][0][2]=1.0*ddsdde[4][0];
  ddsdde4[0][0][2][0]=1.0*ddsdde[4][0];
  ddsdde4[0][0][1][2]=1.0*ddsdde[5][0];
  ddsdde4[0][0][2][1]=1.0*ddsdde[5][0];
  ddsdde4[1][1][0][1]=1.0*ddsdde[3][1];
  ddsdde4[1][1][1][0]=1.0*ddsdde[3][1];
  ddsdde4[1][1][0][2]=1.0*ddsdde[4][1];
  ddsdde4[1][1][2][0]=1.0*ddsdde[4][1];
  ddsdde4[1][1][1][2]=1.0*ddsdde[5][1];
  ddsdde4[1][1][2][1]=1.0*ddsdde[5][1];
  ddsdde4[2][2][0][1]=1.0*ddsdde[3][2];
  ddsdde4[2][2][1][0]=1.0*ddsdde[3][2];
  ddsdde4[2][2][0][2]=1.0*ddsdde[4][2];
  ddsdde4[2][2][2][0]=1.0*ddsdde[4][2];
  ddsdde4[2][2][1][2]=1.0*ddsdde[5][2];
  ddsdde4[2][2][2][1]=1.0*ddsdde[5][2];
  ddsdde4[1][2][0][0]=1.0*ddsdde[0][5];
  ddsdde4[1][2][1][1]=1.0*ddsdde[1][5];
  ddsdde4[1][2][2][2]=1.0*ddsdde[2][5];
  ddsdde4[1][2][0][1]=1.0*ddsdde[3][5];
  ddsdde4[1][2][1][0]=1.0*ddsdde[3][5];
  ddsdde4[1][2][0][2]=1.0*ddsdde[4][5];
  ddsdde4[1][2][2][0]=1.0*ddsdde[4][5];
  ddsdde4[1][2][1][2]=1.0*ddsdde[5][5];
  ddsdde4[1][2][2][1]=1.0*ddsdde[5][5];
  ddsdde4[2][0][0][0]=1.0*ddsdde[0][4];
  ddsdde4[2][0][1][1]=1.0*ddsdde[1][4];
  ddsdde4[2][0][2][2]=1.0*ddsdde[2][4];
  ddsdde4[2][0][0][1]=1.0*ddsdde[3][4];
  ddsdde4[2][0][1][0]=1.0*ddsdde[3][4];
  ddsdde4[2][0][0][2]=1.0*ddsdde[4][4];
  ddsdde4[2][0][2][0]=1.0*ddsdde[4][4];
  ddsdde4[2][0][1][2]=1.0*ddsdde[5][4];
  ddsdde4[2][0][2][1]=1.0*ddsdde[5][4];
  ddsdde4[0][1][0][0]=1.0*ddsdde[0][3];
  ddsdde4[0][1][1][1]=1.0*ddsdde[1][3];
  ddsdde4[0][1][2][2]=1.0*ddsdde[2][3];
  ddsdde4[0][1][0][1]=1.0*ddsdde[3][3];
  ddsdde4[0][1][1][0]=1.0*ddsdde[3][3];
  ddsdde4[0][1][0][2]=1.0*ddsdde[4][3];
  ddsdde4[0][1][2][0]=1.0*ddsdde[4][3];
  ddsdde4[0][1][1][2]=1.0*ddsdde[5][3];
  ddsdde4[0][1][2][1]=1.0*ddsdde[5][3];
  ddsdde4[2][1][0][0]=1.0*ddsdde[0][5];
  ddsdde4[2][1][1][1]=1.0*ddsdde[1][5];
  ddsdde4[2][1][2][2]=1.0*ddsdde[2][5];
  ddsdde4[2][1][0][1]=1.0*ddsdde[3][5];
  ddsdde4[2][1][1][0]=1.0*ddsdde[3][5];
  ddsdde4[2][1][0][2]=1.0*ddsdde[4][5];
  ddsdde4[2][1][2][0]=1.0*ddsdde[4][5];
  ddsdde4[2][1][1][2]=1.0*ddsdde[5][5];
  ddsdde4[2][1][2][1]=1.0*ddsdde[5][5];
  ddsdde4[0][2][0][0]=1.0*ddsdde[0][4];
  ddsdde4[0][2][1][1]=1.0*ddsdde[1][4];
  ddsdde4[0][2][2][2]=1.0*ddsdde[2][4];
  ddsdde4[0][2][0][1]=1.0*ddsdde[3][4];
  ddsdde4[0][2][1][0]=1.0*ddsdde[3][4];
  ddsdde4[0][2][0][2]=1.0*ddsdde[4][4];
  ddsdde4[0][2][2][0]=1.0*ddsdde[4][4];
  ddsdde4[0][2][1][2]=1.0*ddsdde[5][4];
  ddsdde4[0][2][2][1]=1.0*ddsdde[5][4];
  ddsdde4[1][0][0][0]=1.0*ddsdde[0][3];
  ddsdde4[1][0][1][1]=1.0*ddsdde[1][3];
  ddsdde4[1][0][2][2]=1.0*ddsdde[2][3];
  ddsdde4[1][0][0][1]=1.0*ddsdde[3][3];
  ddsdde4[1][0][1][0]=1.0*ddsdde[3][3];
  ddsdde4[1][0][0][2]=1.0*ddsdde[4][3];
  ddsdde4[1][0][2][0]=1.0*ddsdde[4][3];
  ddsdde4[1][0][1][2]=1.0*ddsdde[5][3];
  ddsdde4[1][0][2][1]=1.0*ddsdde[5][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++){
        for (l = 0; l < 3; l++){
        dsdf[i][j][k][l]=0.0;
          for (p = 0; p < 3; p++){
            sum1[i][p][k][l]=0.0;
            sum2[i][p][k][l]=0.0;
            sum3[i][j][k][p]=0.0;
            for (q = 0; q < 3; q++){
              sum11[q][p][k][l]=0.0;
              for (m = 0; m < 3; m++){
                sum11[q][p][k][l]+=
                 det*ddsdde4[q][p][k][m]*finv[l][m]+
                 0.5*det*stre[m][p]*Indim[q][k]*finv[l][m]+
                 0.5*det*stre[q][m]*Indim[p][k]*finv[l][m];};
              sum11[q][p][k][l]-=0.5*det*stre[k][p]*finv[l][q];
              sum11[q][p][k][l]-=0.5*det*stre[q][k]*finv[l][p];
              sum2[i][p][k][l]+=det*stre[q][p]*finv[l][q]*finv[i][k];
              sum3[i][j][k][p]+=det*stre[q][p]*finv[i][q]*finv[j][k];
              sum1[i][p][k][l]+=finv[i][q]*sum11[q][p][k][l];};
            dsdf[i][j][k][l]+=sum1[i][p][k][l]*finv[j][p]-
             finv[j][p]*sum2[i][p][k][l]-sum3[i][j][k][p]*finv[l][p];};};};};};
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      Jac[0][3*i+j]=1.0*dsdf[0][0][i][j];
      Jac[1][3*i+j]=1.0*dsdf[1][1][i][j];
      Jac[2][3*i+j]=1.0*dsdf[2][2][i][j];
      Jac[3][3*i+j]=(dsdf[1][2][i][j]+dsdf[2][1][i][j])/2.0;
      Jac[4][3*i+j]=(dsdf[0][2][i][j]+dsdf[2][0][i][j])/2.0;
      Jac[5][3*i+j]=(dsdf[0][1][i][j]+dsdf[1][0][i][j])/2.0;};};
	  
  return 0;}
