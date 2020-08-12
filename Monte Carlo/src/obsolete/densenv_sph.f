c     ***********************************************************
      subroutine densenv(rad,thetin,costin,sina,pi,pihalf,phi
     1     ,dens)
      
c     calculates density in envelope. called during grid setup
c     so doesn't need to be optimized.

c  history:
c  00/09/06 (baw): write subroutine, modified from opacenv.f
c  01/02/17 (baw): combine opacitybub subroutine into here.
c  03/07/18 (baw): spheroidal envelope

      use tts_mod
      use opacin_mod
      implicit none

      real*8 rad,dens,cosa,rado,radi,phi,thet,thetin,pihalf,y
     $     ,pi,r2,costin,m,n,fact,sina,xmu,rx,rp,zup,xmu0,factor
     $     ,rx2,xmu0new,zlo,zp,f

      integer iflag

c ... begin g77
      real*8 cosd,sind,tand,a
      external cosd,sind,tand
c     ... end g77

c      f=0.d0
c      a=0
c      f=31.60696
c      dens=rhoconst1*(rad/rddust)**a/(1+f**2*costin**2)
      dens=rhoamb
      if (rad.lt.rmine) dens=0.d0
      if (rad.gt.rmax) dens=0.d0

c     print*,'rad,rad/rddust',rad,rad/rddust

      return 
      end

