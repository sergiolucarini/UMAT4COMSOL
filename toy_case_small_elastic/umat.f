      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      include 'aba_param.inc'
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
      integer ii,jj
      real EE,nu,llambda,mu
      EE=props(1)
      nu=props(2)
      llambda=EE*nu/((1.0d0+nu)*(1.0d0-2.0d0*nu))
      mu=EE/(2.0d0*(1.0d0+nu))
      do ii=1,ntens
        do jj=1,ntens
          ddsdde(ii,jj)=0.0d0
        enddo
      enddo
      do ii=1,ndi
        do jj=1,ndi
          ddsdde(ii,jj)=llambda
        enddo
        ddsdde(ii,ii)=llambda+2.0d0*mu
      enddo
      do ii=ndi+1,ntens
        ddsdde(ii,ii)=mu
      enddo 
      do ii=1,ntens
c        stress(ii)=0.0d0
        do jj=1,ntens
          stress(ii)=stress(ii)+ddsdde(ii,jj)*
     1     (dstran(jj))
c     1     (stran(jj)+dstran(jj))
        enddo
      enddo
      return
      end subroutine

