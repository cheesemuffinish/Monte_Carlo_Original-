      subroutine disk

c     monte carlo radiative transfer for disk surrounding star.
c     uses cartesian coordinates, arbitrary density distribution,
c     arbitrary disk structure.
c     This program works for extremely (geometrically) thin disks of large
c     radial extent plus tenuous envelope.  In order to do the
c     radiative transfer in the thin disk, the opacity is calculated
c     at each step, and the photon path integration is variable.
c     In the disk, the step size is small for directions perpendicular
c     to the disk, say zmindisk/10.  Outside the
c     disk, z gt. zmaxdisk, the step size is much larger, zmax/100
c     or zmax/200, something on that order.
c     Note that zmindisk is probably about .1 stellar radius.
c     zmax is about 10000 stellar radii, that is, about 100 AU.
c     The envelope extend is even larger for protostars--10^4 AU,
c     so the step size increases with distance from the source.

c     calls these subroutines:
c        stokes
c        initp
c        opacdisk
c        opacinfall
c        dels
c
c history:
c 99/04/12 (mjw): add cputime for photon loop update output
c                 (this breaks the overall CPU clock in newtts)
c 99/11/20 (baw):  big changes...
c
c     *************************************************************

      implicit none

      include 'stokes.txt'
      include 'tts.txt'
      include 'opacin.txt'
      include 'taunum.txt'
      include 'random.txt'
      include 'out.txt'
      include 'spot.txt'
      include 'vger.txt'
      include 'grid.txt'
      include 'dust.com'
      include 'filt.txt'

      real*8 sini(nmu),xran,t,nu,bnu,zmaxi,xmax,angle,dsl,rnri,rnzi
     1   ,fractx,area,rstarcgs,dconst,rin,rout,dr,thetb
     1   ,zout,zin,dz,rad,normstar1,normstar2,tsum,coslp,sinlp
     1   ,cosnorm,sipnew,fact,rtmp,taunoam,eps
     1   ,rmincgs,tcount,tau_therm,accfrac1,tplanck

      integer ii(3),ns,nd,ne,icount,ithet,gridflag,iphot,i,nstart,ir
     1   ,nend,iz,outflag,it,ip,npacc,ips,idust2,id,npaccstar,iplanck
      logical ifound
      real ran2

c cpu time variables
      real etime,cpusec,time(2)
      character cmsgnm*70
      external etime
     
      write(6,*) 'check, pi',pi
      icount=0

c     set up density grid
      call gridset

c     set up filter functions if making peeled images
      if (ipeel.eq.1) call filt

c     The flux is summed to image(ix,iy,it)
c     polarization to imagei(ixp,iyp,it),imageq...,imageu...

c     inclination arrays.
      do ithet=1,nmu
         sini(ithet)=sqrt(1.d0-u(ithet)**2)
         print*,sini(ithet),nmu
      end do

c     if the inner hole is large, star is effectively point source
c     and theta,phi grid cells are really close together at rmin,
c     so emit photons into radial direction so grid can handle it.
      if (rddust.gt.100) then
         print*,'emitting from star as point source' 
         ips=1
      else
         ips=0
      end if

c     if (ispot.eq.1.and.spotflag.eq.1) call spotset(npsav,wave,tstar)
c     only call spotset once, it will set spotflag=0 after first call

      nscat=0.d0
      tot=0.d0

c     origin of photons, set to 0 for now...
      ior=0

c     x,y arrays (images)
      fractx=2.d0*rmaxi/dble(nx)
c     x,y arrays for peeled-off image
      fractxh=2.d0*rmaxi/dble(nxhst)
      xmax=rmaxi

      rmincgs=rstar*rsol

c     eps for stepping through grid (to step just beyond cell)
c     this number is optimized for the most opaque disk this code
c     can handle right now.
      eps=5.d-8

      call flush(6)
      
c     np=npsav
      print*,'np',np
      flux=0.d0
      sflux=0.d0
      aflux=0.d0
      abflux=0.d0
      scount=0.d0
      dscount=0.d0

c     if (itherm.eq.0) then
c      ns=np+npout
c      nd=0
c      ne=0
c     else
c     calculate luminosity of star
c      nu=2.99792458d14/wave     !in microns
c      call plancknu(tstar,nu,bnu)
c      print*,'stellar Bnu ',bnu
c     bnu=bnu*4*pi     !convert from energy/s/cm2/Sr to energy/s/cm2
c      rstarcgs=rstar*rsol
c      area=pi*(rstarcgs)**2
c      normstar=area*bnu         !this is an energy/s; normalization 
c                                !for output
c      print*,'stellar energy at this wavelength at earth',normstar
c      print*,'luminosity ',4.d0*area/3.826d33*tstar**4*5.669d-5 
c                                !energy/s
c      normstar1=normstar*4.d0*pi !flux=pi*bnu area=4*pi*rstar**2
c      print*,'luminosity at this wavelength ',normstar1
      ns=(np+npout)*(1.d0-accfrac)
      npacc=(np+npout)*accfrac
      np=ns+npacc+npout
      npaccstar=accfrac2*ns
      print*,'ns,np',ns,np
      print*,'np_disk,np_shock',npacc,npaccstar

c      ns=np+npout
c     end if

c     convert Tstar to Jon's units
      tstar=tstar/11605.d0
      tshock=tshock/11605.d0

c     ********  do loop over stellar photons  *******
      print*,'ns',ns
      print*,'npout',npout

      do iphot=1,np
         
c     if (iphot.gt.30000.and.iphot.lt.31000) then
c     continue
c     print*,iphot
c     end if
         if(mod(iphot,iwrite).eq.0) then
c     cpu time
            cpusec=etime(time)
c     use explicit format to keep g77 happy
            write(cmsgnm,'(i12,a20,f11.2,a)') iphot,
     $           ' photons completed. ',cpusec,' = CPU time (sec)'
            call WRIMSG('MAIN',cmsgnm)
         end if
         
c         if (ran().lt.accfrac) then
         if (iphot.lt.(npacc)) then
c     disk accretion photons
            call initpacc(ii,eps,ns)
            call opacset(nub)
c            print*,'rsq.rtot,zp',rsq,rtot,zp,ii(1),ii(2),ii(3)
            wave=1.2398d0/nub
c     wave=2.99792458d14/nub
c     print*,'nub',nub
c     print*,'wave',nub
	      ir=ii(1)
	      it=ii(2)
	      ip=ii(3)
	      idust2=dustarr(ir,it,ip)
	      do id=1,4
                 kapd(id)=kappa(id)*rmincgs*kappav(id)
                 rlam(id)=albedo(id)
              end do
c     print*,'kappa,kappav,kapd',kappa,kappav,kapd
            iflag=0
            if (idust2.eq.1.or.idust2.eq.2) then
               ior=2
            else
               ior=3
            end if
c            ior=2
            if (ipeel.eq.1) then
               call peeloff(xp,yp,zp,sipnew,sqp,sup,svp
     1              ,cost,sint,cosp,sinp,phi
     1              ,pi,r2p,hit,htot,rsq,rtot
     1              ,tsum,ii,idust,iflag,iphot,eps)
            end if
c******  check peeling off of accretion photon!!!!!!!  071902
c        set normalization!!!!!

         else if (iphot.lt.(npacc+ns)) then
c     star illumination
            if (iphot.lt.(npacc+npaccstar)) then
               iplanck=1
               if (ips.eq.0) then
                  call initp(iplanck)
               else
                  call initp_ps(iplanck)
               end if
c     emit half the photons as x-rays
               if (ran().lt.0.5) then
                  wave=0.01+ran()*0.05
                  nub=1.2398d0/wave
               end if
            else
               iplanck=0
               if (ips.eq.0) then
                  call initp(iplanck)
               else
                  call initp_ps(iplanck)
               end if
            end if
c     origin of photon ior = 1 for star
            ior=1
            call opacset(nub)
            wave=1.2398d0/nub
c     print*,'nub',nub
c     print*,'wave',nub
            idust2=dustarr(1,1,1)
            do id=1,4
               kapd(id)=kappa(id)*rmincgs*kappav(id)
               rlam(id)=albedo(id)
            end do
c     print*,'kappa,kappav,kapd',kappa,kappav,kapd
c     print*,'cost,sint',cost,sint
            coslp=cos(lp)
            sinlp=sin(lp)
c     zstar
            zp=cosb*rmin*(1.d0+eps)
c     rstar=xp if all photons start off on x-axis
            rp=sinb*rmin*(1.d0+eps)
            xp=rp*coslp
            yp=rp*sinlp
            rsq=zp**2+rp**2
            rtot=sqrt(rsq)
            ux=sint*cosp
            uy=sint*sinp
            uz=cost
            iflag=0
            ii(1)=1             !index of radial grid (stellar 
                                !surface is 1)
c     if (ntg.gt.1) then
            thetb=acos(cosb)
            call locate(thetarr,ntg,thetb,it)
            ii(2)=it
c     else
c     ii(2)=1
c     it=1
c     end if
c     if (npg.gt.1) then
            call locate(phiarr,npg,lp,ip)
            ii(3)=ip
c     else
c     ii(3)=1
c     end if
            
c     print*,'xp,yp,zp,ux,uy,uz',xp,yp,zp,ux,uy,uz
c     print*,'ir,it,ip',ii(1),ii(2),ii(3)
c     print*,'cost,sint,cosp,sinp',cost,sint,cosp,sinp
c     stop
c     first, peel off direct flux
c     weight photon intensity by angle between normal and 
c     photon direction, since emitted from a surface.
            if (ipeel.eq.1) then
                  cosnorm=cosb*coste+(sinb*sinte*(cospe*coslp+
     1                 sinpe*sinlp))
                  if (limb.eq.0) then
c     intensity constant, energy/Sr proportional to mu
                     sipnew=4.*sip*cosnorm !normalization = 2
                  else
c     intensity goes as (1+mu).  energy/Sr has another factor
c     of mu.
                     sipnew=12.d0/5.d0*sip*(cosnorm+cosnorm*cosnorm)
                  end if
c     sipnew=sip	 
                  if (cosnorm.gt.0.d0) then
                     call peeloff(xp,yp,zp,sipnew,sqp,sup,svp
     1                    ,cost,sint,cosp,sinp,phi
     1                    ,pi,r2p,hit,htot,rsq,rtot
     1                    ,tsum,ii,idust,iflag,iphot,eps)   
                  end if
c     vger stuff is not properly implemented to radeq code
c     but leaving it in
c     for future use when it is implemented.
c     for now, iveeg=0
               if (cosnorm.gt.0.and.iveeg.eq.1) then
                  call vger(xp,yp,zp,sipnew,sqp,sup,svp
     1                 ,cost,sint,cosp,sinp,phi
     1                 ,pi,r2p,hit,htot,rsq,rtot
     1                 ,tsum,ii,idust,iflag,iphot,eps)   
               end if
            end if
         else
c     outside illumination
            call initpout
            ux=sint*cosp
            uy=sint*sinp
            uz=cost
            iflag=0
c     don't peel off direct flux because star is outside image
c     field  . 
c     but send star flux to vger   
c     vger stuff is not properly implemented to radeq code
c     but leaving it in
c     for future use when it is implemented.
c     for now, iveeg=0
            if (iveeg.eq.1) then
               sipnew=sip
               call vger(xse,yse,zse,sipnew,sqp,sup,svp
     1              ,cost,sint,cosp,sinp,phi
     1              ,pi,r2p,hit,htot,rsq,rtot
     1              ,tsum,ii,idust,iflag,iphot,eps)   
               if (iphot.eq.1) then
                  print*, 'tau from outside star to vger ',
     1                 tsum
                  write(12,*) 'tau from out star to vger ',
     1                 tsum
               end if
            end if
c     now calculate position in envelope of photon of selected
c     direction. note that it may not hit envelope at all.
            call Rdist(ifound,t,ux,uy,cost,xs,ys,zs,rmax)
            if (ifound.eqv..false.) go to 5
            fact=1.0d0+1.d-3
            xp=xs+ux*t*fact
            yp=ys+uy*t*fact
            zp=zs+cost*t*fact
            rsq=xp**2+yp**2+zp**2
            rtot=sqrt(rsq)
c     rp=sqrt(xp**2+yp**2)
c     rtmp=sqrt(xp**2+yp**2+zp**2)	 
         end if
c     now calculate scattered flux
c     print*,'before propagate, ir,it,ip',ii(1),ii(2),ii(3)
         call propagate(iphot,sini,fractx,ii,eps,ns,idust2
     $        ,icount)
c     print*,'done with propagate'
   
 5    end do
      
c     ***** end of loop over stellar photons  ******
      
c      np=ns

c     calculate final temp using lucy method
      call tfinal(ns)

      write(6,*) 'sini,cosi'
      write(6,*) (sini(i),i=1,nmu),(u(i),i=1,nmu)
      write(6,*) 'fractx,xmax',fractx,xmax
      
      write(6,*) 'fraction of phots scattered outside image'
      write(6,*) dble(icount)/dble(np)
      print*,np,scount,dscount
      tcount=(dble(np)-scount-dscount)/dble(np)
      tau_therm=-log(1.d0-tcount)
      write(12,*) 'tcount, tau_therm',tcount, tau_therm
      
      return
      end

c     *********************************************************
