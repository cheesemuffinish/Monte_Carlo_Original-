      real kappa(4),kappav(4),albedo(4),g(4),pl(4),pc,sc,g2(4)
     $	,onemg2(4),g2p1(4)
     +     ,twog(4),p1maxinv(4),Tcond
      common /dustparm/kappa,kappav,albedo,g,pl,pc,sc,g2,
     +                 onemg2,g2p1,twog,p1maxinv,Tcond
      
      integer nlammax,nlambda(4)
      parameter (nlammax=200)
      real lamdust,nudust,kappad,adust,gdust,pldust      
      dimension lamdust(nlammax,4),nudust(nlammax,4),kappad(nlammax,4),
     +          adust(nlammax,4),gdust(nlammax,4),pldust(nlammax,4)
      common /dustarrays/lamdust,nudust,kappad,adust,gdust,pldust,
     +                   nlambda

      integer nT,nnu
      parameter(nT=1000,nnu=1025)
      real lTint,Tint,nuint,kappap,kappar,kapBint,kapdBint,dBint,kappaf
      dimension lTint(nT),Tint(nT),nuint(nnu),kappap(nT,4),kappar(nT,4),
     +          kapBint(nT,nnu,4),kapdBint(nT,nnu,4),dBint(nT,nnu,4),
     +	    kappaf(4)
      common /kB_interp/lTint,Tint,nuint,kappap,kappar,
     +                  kapBint,kapdBint,dBint,kappaf
