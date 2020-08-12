      real tauv,tau0,nsigma,rmax,r2max,h0,h02,mu0,Asphere,alpha,
     +     alphap,beta,betam1,bp2ma,mdisk,sigmu0,tauFpole
      common /densparm/tauv,tau0,nsigma,rmax,r2max,h0,h02,mu0,
     +                 Asphere,alpha,alphap,beta,betam1,bp2ma,
     +                 mdisk,sigmu0,tauFpole

      integer nr,nr1,nmu,nmu1,nc
      parameter (nr=100,nr1=nr+1,nmu=30,nmu1=nmu+1,nc=nr*nmu)
      integer nrgrid,nabs
      real Nrfac,rbin,mubin,taumid,mass,Tdust,Tdisk
      dimension rbin(nr1),taumid(nr),mubin(nr,nmu1),mass(nr,nmu),
     +          Tdust(nr,nmu),nabs(nr,nmu),Tdisk(nr)
      common /dustgrid/Nrfac,rbin,taumid,mubin,mass,Tdust,Tdisk,
     +                 nrgrid,nabs
