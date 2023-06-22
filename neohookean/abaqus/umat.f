C     UMAT Neo-Hookean hyperelastic model
c     umat.f
c     Copyright (C) 2023  Sergio Lucarini
c
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      real*8 bbar(6), distgr(3,3),cE,cnu,c10,d1,det
      real*8 cG,cK,pr,trbbar,cG23
      integer ii,jj

      cE=PROPS(1)
      cnu=PROPS(2)
      c10 = cE / (4. * (1. + cnu))
      d1 = 6. * (1. - 2. * cnu) / cE
      det=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)
     1   -DFGRD1(1, 2)*DFGRD1(2, 1)*DFGRD1(3, 3)
     2   +DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)
     3   +DFGRD1(1, 3)*DFGRD1(3, 2)*DFGRD1(2, 1)
     4   -DFGRD1(1, 3)*DFGRD1(3, 1)*DFGRD1(2, 2)
     5   -DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
      do ii=1, 3
         do jj=1, 3
           distgr(jj, ii)=det**(-1.0d0/3.0d0)*DFGRD1(jj, ii)
         enddo
      enddo
      bbar(1)=distgr(1, 1)**2+distgr(1, 2)**2+distgr(1, 3)**2
      bbar(2)=distgr(2, 1)**2+distgr(2, 2)**2+distgr(2, 3)**2
      bbar(3)=distgr(3, 3)**2+distgr(3, 1)**2+distgr(3, 2)**2
      bbar(4)=distgr(1, 1)*distgr(2, 1)+distgr(1, 2)*distgr(2, 2)
     1       +distgr(1, 3)*distgr(2, 3)
      bbar(5)=distgr(1, 1)*distgr(3, 1)+distgr(1, 2)*distgr(3, 2)
     1         +distgr(1, 3)*distgr(3, 3)
      bbar(6)=distgr(2, 1)*distgr(3, 1)+distgr(2, 2)*distgr(3, 2)
     1         +distgr(2, 3)*distgr(3, 3)
      trbbar=(bbar(1)+bbar(2)+bbar(3))/3.0d0
      cG=2.0d0*c10/det
      cK=2.0d0/d1*(2.0d0*det-1.0d0)
      pr=2.0d0/d1*(det-1.0d0)
      do ii=1,3
         STRESS(ii)=cG*(bbar(ii)-trbbar)+pr
      enddo
      do ii=4,6
         STRESS(ii)=cG*bbar(ii)
      enddo
      cG23=cG*2.0d0/3.0d0
      DDSDDE(1, 1)= cG23*(bbar(1)+trbbar)+cK
      DDSDDE(2, 2)= cG23*(bbar(2)+trbbar)+cK
      DDSDDE(3, 3)= cG23*(bbar(3)+trbbar)+cK
      DDSDDE(1, 2)=-cG23*(bbar(1)+bbar(2)-trbbar)+cK
      DDSDDE(1, 3)=-cG23*(bbar(1)+bbar(3)-trbbar)+cK
      DDSDDE(2, 3)=-cG23*(bbar(2)+bbar(3)-trbbar)+cK
      DDSDDE(1, 4)= cG23*bbar(4)/2.0d0
      DDSDDE(2, 4)= cG23*bbar(4)/2.0d0
      DDSDDE(3, 4)=-cG23*bbar(4)
      DDSDDE(4, 4)= cG*(bbar(1)+bbar(2))/2.0d0
      DDSDDE(1, 5)= cG23*bbar(5)/2.0d0
      DDSDDE(2, 5)=-cG23*bbar(5)
      DDSDDE(3, 5)= cG23*bbar(5)/2.0d0
      DDSDDE(1, 6)=-cG23*bbar(6)
      DDSDDE(2, 6)= cG23*bbar(6)/2.0d0
      DDSDDE(3, 6)= cG23*bbar(6)/2.0d0
      DDSDDE(5, 5)= cG*(bbar(1)+bbar(3))/2.0d0
      DDSDDE(6, 6)= cG*(bbar(2)+bbar(3))/2.0d0
      DDSDDE(4,5)= cG*bbar(6)/2.0d0
      DDSDDE(4,6)= cG*bbar(5)/2.0d0
      DDSDDE(5,6)= cG*bbar(4)/2.0d0
      do ii=1,6
         do jj=1,ii-1
            DDSDDE(ii, jj)=DDSDDE(jj, ii)
         enddo
      enddo
      RETURN
      END
