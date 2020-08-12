      subroutine gridset

c  set up grid:  3-D spherical grid with variable spacing in r and 
c  theta (set by exponents, rexp, texp)
c  make a bunch of arrays necessary for find_wall subroutine

c history:
c 00/09/06 (baw): write subroutine
	
      use tts_mod
      use grid_mod
      use opacin_mod
      use out_mod
      implicit none
      use dust_mod

      real*8 dr,dt,dp,rad,phi,thet,pi,r2p,densd,dense,vol
     $ ,cost,sint,pihalf,eps,tau,rming1,rming2,rming3,tiny
     $     ,rmincgs,taud,taue,mu0,tauv,taufmin,mu,thetatm
     $     ,thetmin,tauzrat,sigmu,rsonr,Tdisk,taudisk
     $     ,mudisk,c1,a1,tauw,rtop,rfac,tauRr,taumid,tauRth
     $     ,tauRos,res,thettop,maxdens,tauRmu,radamb,rbeg2
      real xyarr(200,200),x,y,xmax,dx,r,dens1
      real chiR,ierfc,erfc
      integer ir,it,ip,nttmp,gflag,n,ix,iy,nr0,ntatm,nrwall,irbeg,
     $     dcount,id,ide,idd
      external erfc

      pi=4.d0*datan(1.d0)
      r2p=2.d0*pi
      pihalf=pi/2.d0
      tiny=1.d-15
      rmincgs=rstar*rsol
      radamb=rmax

c     zero arrays
      do ir=1,nrg
      do it=1,ntg
      do ip=1,npg
         densarr(ir,it,ip)=0.d0
         massarr(ir,it,ip)=0.d0
      end do
      end do
      end do
c     convert opacity to units of cm^2/gm*rstar(cgs) because distance
c     is in units if 1/rstar, and dtau=kapd*rho*ds
      do id=1,ndg
         kapd(id)=kappav(id)*rmincgs
      end do
c     print*,'kappav',kappav
c     stop

c     calculate some constants for the TSC envelope density
      call envset

c     make grid include minimum and maximum values.

      print*, ' '
      print*, 'grid setup'
      pihalf=pi/2.d0

      print*,'rddust,rmine,rmin',rddust,rmine,rmin
c     rgrid
      gflag=0
      rarr(1)=rmin
      if (rddust.lt.rmine) then
         gflag=1
         rarr(2)=rddust
         rming1=rddust
         rming2=rmine
         irbeg=2
      else if (rddust.gt.rmine) then
         gflag=2
         rarr(2)=rmine
         rming1=rmine
         rming2=rddust
         irbeg=3
      else if (rddust.eq.rmine) then
         gflag=3
         rarr(2)=rmine
         rming1=rmine
         irbeg=2
      end if
      if (massd.eq.0.d0) then
c     since there's empty space below rmine and rddust....
         dr=(rmax-rming1)/(dble(nrg-2))**rexp
         print*,'rexp,dr,rmin',rexp,dr/autors,rmin/autors
         do ir=3,nrg
            rarr(ir)=rming1+dr*(dble(ir-2))**rexp
c     print*,'rarr ',ir,rarr(ir)/au
         end do
         
      else
c     lots of fine spacing in inner region of disk
         
         a1=1.d0-a
         
         mu0=z1*Rddust**b/Rddust
         tauv=taur
         tauzrat=sqrt(pihalf)*a1*mu0/((rmaxd/rddust)**a1-1.d0)
         tauRos=tauv*chiR(1600.d0/11605.d0,1)   !setting idust2=1, disk
         tauRth=tauzrat*tauRos
c         taudisk=10.0
         taudisk=min(10.d0,tauRth)
         
         print*,'tauRos,tauv,tauRth',tauRos,tauv,tauRth
         tauw=0.5
         dr=(rmax-rming1)/(dble(nrg-irbeg))**rexp
         rbeg2=rming1+dr
         print*,'rbeg2',rbeg2
         if (tauRos.gt.0.51) then
         rfac=(1.d0-(tauw*(1.d0-(rmaxd/rddust)**a1))/tauRos)**(1./a1)
         if (rfac.lt.1.0000001d0) then
            print*,''
            print*,'WARNING: not enough precision to resolve inner disk'
            print*,'rfac = ',rfac
c            rfac=rfac*2.d0
            rfac=1.000001d0
             print*,'resetting rfac = ',rfac
            print*,''
         end if

         mu=sqrt(2.0)*ierfc(sngl(taudisk/tauRth))
         rtop=Rddust*(1.d0-(taudisk*(1.d0-(rmaxd/rddust)**a1))/
     +        (tauRos*exp(-0.5d0*mu*mu))      )**(1./a1)
c     +        (tauRos)      )**(1./a1)
         
         print*,'Rddust,taudisk,rmaxd/rddust,tauRos,a1',
     1        Rddust,taudisk,rmaxd/rddust,tauRos,a1
         print*,'rfac,rtop,mu',rfac,rtop,mu
         ir=irbeg
         rarr(ir)=rddust
         if (rfac*rddust.lt.rbeg2) then
         print*,'making gaussian spacing in r, rtop',rtop
         do while ((rarr(ir).lt.rtop).and.(ir.lt.(nrg-1)))
            rarr(ir+1)=rfac*rarr(ir)
            rfac=rfac**1.1
            ir=ir+1
c            print*,'rarr',rarr(ir),rtop,ir,rfac
         end do
         if(ir.ge.nrg) then
            write(*,*) 'too many points in rgrid'
            stop
         end if
         end if
         else
	   ir=irbeg
         end if

         nrwall=ir
         print*,'nrwall,rarr(nrwall)',nrwall,rarr(nrwall)
         
         rming1=rarr(nrwall)
         dr=(rmax-rming1)/(dble(nrg-nrwall))**rexp
         print*,'rexp,dr,rmin',rexp,dr/autors,rming1
         do ir=nrwall+1,nrg
            rarr(ir)=rming1+dr*(dble(ir-nrwall))**rexp
         end do
         
      end if

c     set rmine to a grid location
      if (gflag.eq.1) then
         call locate(rarr,nrg,rmine,ir)
         rarr(ir+1)=rmine
      end if

c     r-squared array
      do ir=1,nrg
         r2arr(ir)=rarr(ir)**2
c     print first few points of rarr
         if (ir.lt.10) then
            print*,'rarr/rmin,ir ',rarr(ir)/rmin,ir
         end if
      end do

      open(unit=15,file='rarr.dat',status='unknown')
      write(15,*) nrg,' = number of grid points in r'
      write(15,*) 'index     r/rstar  r(au)'
      do ir=1,nrg
         write(15,*) ir,(rarr(ir)),sngl(rarr(ir)/autors)
      end do
      close(15)

      print*,'rarr(nrg),r2arr(nrg)',rarr(nrg),r2arr(nrg)

      if (ntg.gt.1) then
c     set up a tmp array going from 0 - 90 from equ. to pole
c     make theta=0 bin 5 degrees wide (otherwise, too much noise,
c     nothing happens there anyway)
         thettop=5.*deg2rad
         if (massd.eq.0.d0) then

            ntatm=1
            thetatm=0.d0
            tmptharr(1)=0.d0
            nttmp=(ntg+1)/2

         else

c     taufmin=min(0.1,0.1*kappaf*tauv)
c	      use idust2=2, for kappaf, for disk dust properties
            taufmin=min(.001d0,0.001d0*kappaf(2)*tauv)
            
            print*,'kappaf*tauv',kappaf(2)*tauv
            
            mu=min(mu0*sqrt(2.d0*log(tauv*kappaf(2)/taufmin)),0.5d0)
            print*,'tauv,kappaf,taufmin',tauv,kappaf(2),taufmin
            

c     ************
c     for disk-only model, let mu=0.5 to sample entire disk height
c     at high-res
c     mu=0.5
c     ***********

            thetatm=pihalf-acos(mu)
c     test
            print*,'thetatm',thetatm*rad2deg
            nttmp=(ntg+1)/2

c     *********************
c     change this for disk or envelope runs
c     for envelope with 1 degree resolution in polar region
c     except for first bin (theta=5).
            res=1.
c     for disk with nothing in the polar region
c     res = size of angle bin in poles, approximately 
c     res=10.
c     *********************

            ntatm=nttmp-(pihalf-thetatm-thettop)*90./pihalf/res

            thetmin=min(0.1d0*asin(mu0),thetatm/ntatm)
            
            tmptharr(1)=0.d0
            do it=2,ntatm
               tmptharr(it)=thetmin *
     1              (thetatm/thetmin)**((it-2.d0)/(ntatm-2.d0))
c     print*,'tmptharr ',tmptharr(it)*rad2deg
            end do            
         end if
         do it=ntatm+1,nttmp-1
            tmptharr(it)=thetatm+
     1           (it-ntatm)*(pihalf-thetatm-thettop)/(nttmp-ntatm)
c            print*,'tmptharr(it)',tmptharr(it)
         end do
         tmptharr(nttmp)=pihalf
         
         thetarr(nttmp)=pihalf
         do it=2,nttmp-1
            thetarr(nttmp+1-it)=pihalf-tmptharr(it)
            thetarr(nttmp-1+it)=pihalf+tmptharr(it)
         end do
         thetarr(1)=0.d0
         thetarr(ntg)=pi
         
c     this is not the eps used in the rest of the code!  see
c     newdisktherm for that
         eps=1.d-8
         do it=1,ntg
            if (it.eq.1) then
               thetarr(it)=0.d0
               costarr(it)=1.d0
               sintarr(it)=0.d0
               tan2arr(it)=0.d0
            else if (it.eq.ntg) then
               thetarr(it)=pi
               costarr(it)=-1.d0
               sintarr(it)=0.d0
               tan2arr(it)=0.d0
            else if (it.eq.nttmp) then
               thetarr(it)=pihalf
               costarr(it)=0.d0
               sintarr(it)=1.d0
               tan2arr(it)=-1.d0
            else
               costarr(it)=cos(thetarr(it))
               sintarr(it)=sin(thetarr(it))
               tan2arr(it)=tan(thetarr(it))**2
            end if
c            print*,'thetarr,costarr,tan2arr '
c            print*,'thetarr',thetarr(it)*rad2deg
         end do
      else
         thetarr(1)=0.d0
         costarr(1)=1.d0
         sintarr(1)=0.d0
         tan2arr(1)=0.d0
      end if

      open(unit=15,file='tharr.dat',status='unknown')
      write(15,*) ntg,' = number of grid points in theta'
      write(15,*) 
     1     'index  thet(rad)  thet(deg)    cost     tan**2(thet)'
      do it=1,ntg
         write(15,*) it,sngl(thetarr(it)),sngl(thetarr(it)*rad2deg),
     1        sngl(costarr(it)),sngl(tan2arr(it))
      end do
      close(15)

      if (npg.gt.1) then
         dp=r2p/dble(npg-1)
         do ip=1,npg
            phiarr(ip)=dp*(dble(ip-1))
            aarr(ip)=sin(phiarr(ip))
            barr(ip)=-cos(phiarr(ip))
c     print*,'phiarr, aarr, barr ',phiarr(ip),aarr(ip),barr(ip)
c     carr(ip)=
c     darr(ip)=
         end do
      else
         dp=r2p
         ip=1
         phiarr(ip)=0.d0
         aarr(ip)=sin(phiarr(ip))
         barr(ip)=-cos(phiarr(ip))
      end if

      open(unit=15,file='phiarr.dat',status='unknown')
      write(15,*) npg
      do ip=1,npg
         write(15,*) phiarr(ip)*rad2deg
      end do
      close(15)

c     set up diffusion grid
      dcount=0
      if (massd.eq.0.d0) then
         do ir =1,nrg
         do it=1,ntg
         do ip=1,npg
            diffus(ir,it,ip)=.false.
            diffdir(ir,it,ip)=0
         end do
         end do
         end do
      else
         if (rddust.eq.rmin) then 
            nr0=1
         else
            nr0=3
            do it=1,ntg
               do ip=1,npg
                  diffus(1,it,ip) = .false.
                  diffdir(ir,it,ip)=0
               end do
            end do
         end if
         do ir=nr0,nrg-1
            r=0.5d0*(rarr(ir)+rarr(ir+1))/rddust
            sigmu=mu0*r**(b-1.d0)
            Rsonr=1./(rddust*r)
            Tdisk=min(1600.d0/11605.d0,max(3.d0/11605.d0,(Tstar/11605.d0)
     +           *(max(2./3.*(Rsonr)**3,
     +           (asin(Rsonr)-(Rsonr)*sqrt(1.d0-(Rsonr)**2)))/pi
c     +        +   (3.*Lacc*(Rsonr)**3*(1.d0-sqrt(Rsonr)))
     +           )**0.25))
            taumid=tauzrat*tauv*chiR(sngl(Tdisk),1)/r**(a-b)
            if (taumid.gt.taudisk) then
               mudisk=sigmu*sqrt(2.d0)*ierfc(sngl(taudisk/taumid))
            else
               mudisk=0.
            end if
c     print*,'ir,taudisk,taumid,mudisk',ir,taudisk,taumid,mudisk
            if (ntg.gt.1) then
               do it=1,ntg-1
                  mu=abs(cos(0.5d0*(thetarr(it)+thetarr(it+1))))
                  tauRmu=taumid*erfc(sngl(mu/(sigmu*sqrt(2.0))))
                  tauRr=tauRos*exp(-0.5d0*mu*mu/(sigmu*sigmu))*
     +                 (1.d0-r**a1)/(1.d0-(rmaxd/rddust)**a1)
                  do ip=1,npg
                     if ((mu.lt.mudisk).and.(tauRr.gt.taudisk)) then
                        diffus(ir,it,ip) = .true.
                        if (tauRr.lt.tauRmu) then
                           diffdir(ir,it,ip) = -1
                        else
                           if (thetarr(it).lt.pihalf) then
                              diffdir(ir,it,ip) = -2
                           else
                              diffdir(ir,it,ip) = 2
                           end if
                        end if
                        dcount=dcount+1
c     print*,'true'
                     else
                        diffus(ir,it,ip) = .false.
                        diffdir(ir,it,ip) = 0
c     print*,'mu,mudisk,tauRr,taudisk'
c     print*,mu,mudisk,tauRr,taudisk
c     print*,'false'
                     end if
                  end do
               end do
            else
               do it=1,ntg
                  do ip=1,npg
                     diffus(ir,it,ip) = .false.
                     diffdir(ir,it,ip) = 0
                  end do
               end do
            end if
         end do
      end if

c     testing
c      do ir =1,nrg
c         do it=1,ntg
c            do ip=1,npg
c               diffus(ir,it,ip)=.false.
c            end do
c         end do
c      end do
c      dcount=0
c      print*,
c     $ 'WARNING, TURNED OFF DIFFUSION, see line 378-387 of gridset.f'

ctest
c      do ir =1,nrg
c         do it=1,ntg
c            do ip=1,npg
c               if (abs(diffdir(ir,it,ip)).eq.1) then
c                  diffus(ir,it,ip)=.false.
c                  diffdir(ir,it,ip)=0
c               end if
c            end do
c         end do
c      end do 
     
      print*,'number of diffusion cells in grid',dcount

      open(unit=15,file='diffuse.unf',status='unknown',
     1     form='unformatted')
      write(15) nrg,ntg,npg
      write(15) (rarr(ir)/autors,ir=1,nrg)
      write(15) (thetarr(it),it=1,ntg)
      write(15) (phiarr(ip),ip=1,npg)
      write(15) (((diffus(ir,it,ip),ir=1,nrg),it=1,ntg)
     1     ,ip=1,npg)
      close(15)

      open(unit=15,file='diffdir.unf',status='unknown',
     1     form='unformatted')
      write(15) nrg,ntg,npg
      write(15) (rarr(ir)/autors,ir=1,nrg)
      write(15) (thetarr(it),it=1,ntg)
      write(15) (phiarr(ip),ip=1,npg)
      write(15) (((diffdir(ir,it,ip),ir=1,nrg),it=1,ntg)
     1     ,ip=1,npg)
      close(15)

c      print*,'massd',massd
c      stop

c     calculate density in grid
      massenv=0.d0
      massdisk=0.d0
      maxdens=0.d0
      do ir=1,nrg-1
         rad=0.5d0*(rarr(ir)+rarr(ir+1))
         dr=rarr(ir+1)-rarr(ir)
         if (ntg.gt.1) then     !2- or 3-D atmosphere
            do it=1,ntg-1
               thet=0.5d0*(thetarr(it)+thetarr(it+1))
               cost=cos(thet)
               sint=sin(thet)
               dt=thetarr(it+1)-thetarr(it)
               if (npg.gt.1) then
                  do ip=1,npg-1   !3-D atmosphere
                     phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
                     dp=phiarr(ip+1)-phiarr(ip)
                     call densenv(rad,thet,cost,sint,phi
     1                    ,dense,radamb,ide)
                     vol=rad**2*sint*dt*dp*dr*rmincgs**3
                     massenv=massenv+dense*vol
                     call densdisk(rad,sint,cost,phi,densd,idd)
                     massdisk=massdisk+densd*vol
                     densarr(ir,it,ip)=dense+densd
                     if (densarr(ir,it,ip).gt.maxdens) 
     $                    maxdens=densarr(ir,it,ip)
                     massarr(ir,it,ip)=(dense+densd)*vol
                     massarr(ir,it,ip)=massarr(ir,it,ip)**0.25
c                     if (idd.eq.0) then
c                        dustarr(ir,it,ip)=ide
c                     else
                        if (dense.gt.densd) then
                           dustarr(ir,it,ip)=ide
                        else
                           dustarr(ir,it,ip)=idd
                        end if
c                     end if
                  end do
               else           !2-D atmosphere
                  ip=1       
                  phi=0.d0
                  dp=r2p
                  call densenv(rad,thet,cost,sint,phi
     1                 ,dense,radamb,ide)
                  vol=rad**2*sint*dt*dp*dr*rmincgs**3
                  massenv=massenv+dense*vol
                  call densdisk(rad,sint,cost,phi,densd,idd)
                  massdisk=massdisk+densd*vol
                  densarr(ir,it,ip)=dense+densd
                  if (densarr(ir,it,ip).gt.maxdens) 
     $                 maxdens=densarr(ir,it,ip)
                  massarr(ir,it,ip)=(dense+densd)*vol
                  massarr(ir,it,ip)=massarr(ir,it,ip)**0.25
c                  if (idd.eq.0) then
c                     dustarr(ir,it,ip)=ide
c                  else
                     if (dense.gt.densd) then
                        dustarr(ir,it,ip)=ide
                     else
                        dustarr(ir,it,ip)=idd
                     end if
c                  end if
               end if
            end do
         else              !1-D atmosphere
            phi=0.d0
            ip=1
            it=1
c            tdust(ir,it,ip)=1000.d0/11605.d0
            thet=60.d0*deg2rad
            cost=cos(thet)
            sint=sin(thet)
            call densenv(rad,thet,cost,sint,phi,dense,
     $		radamb,ide)
            vol=rad**2*dr*2.d0*r2p*rmincgs**3
            massenv=massenv+dense*vol
            call densdisk(rad,sint,cost,phi,densd,idd)
            massdisk=massdisk+densd*vol
            densarr(ir,it,ip)=dense+densd
            if (densarr(ir,it,ip).gt.maxdens) 
     $           maxdens=densarr(ir,it,ip)
            massarr(ir,it,ip)=(dense+densd)*vol
            massarr(ir,it,ip)=massarr(ir,it,ip)**0.25
c            if (idd.eq.0) then
c               dustarr(ir,it,ip)=ide
c            else
              if (dense.gt.densd) then
                 dustarr(ir,it,ip)=ide
              else
                 dustarr(ir,it,ip)=idd
              end if
c            end if
c     print*,'densarr',densarr(ir,it,ip)
         end if
      end do
      massenv=massenv/msun
      massdisk=massdisk/msun
c     important note:  densarr(ir,it,ip) is the density halfway between
c     ir and ir+1, it and it+1, ip and ip+1; it is NOT the density at
c     ir,it,ip.

      print*, 'massenv (in solar units) ',massenv
      print*, 'massdisk from grid, compared to input ',massdisk,
     $     massd
      print*,'max density in disk/envelope',maxdens
      print*,'min radius where rho=rhoamb ',radamb/autors
      print*, ' '

c     write out rarr,thetarr,densarr at phi=0 so you can read them into
c     IDL.

c	check dustarr
	do ir=1,nrg
	do it=1,ntg
	do ip=1,npg
           if (rarr(ir).lt.rmine.and.rarr(ir).lt.rmind) 
     $          dustarr(ir,it,ip)=4
	   if (dustarr(ir,it,ip).eq.0) then
	     if (ir.ne.nrg.and.it.ne.ntg) then
	        print*,'oops, dustarr,r,th,phi',
     $	     rarr(ir),thetarr(it),phiarr(ip)	
            end if
            dustarr(ir,it,ip)=4    !outflow dust    
	   end if
	end do
	end do
	end do

c     integrate optical depth along the theta angles
c     at phi=0.
      open(unit=15,file='tau.dat',status='unknown')
      write(15,*) 
     1 '    theta          tau_env        tau_disk       tau_tot'
      do it=1,ntg-1
         tau=0.d0
         taud=0.d0
         taue=0.d0
         if (ntg.eq.1) then
            thet=60.d0*deg2rad
         else
            thet=0.5d0*(thetarr(it)+thetarr(it+1))
         end if
         if (ntg.gt.1.and.it.eq.(ntg/2)+1) then
c     TSC blows up at thet=90
            print*,'thet',thet*rad2deg
            thet=89.999d0*deg2rad
         end if
         cost=cos(thet)
         sint=sin(thet)
         phi=0.d0
         do ir=1,nrg-1
            dr=rarr(ir+1)-rarr(ir)
            rad=0.5d0*(rarr(ir)+rarr(ir+1))
            call densenv(rad,thet,cost,sint,phi,dense,
     &		radamb,ide)
            call densdisk(rad,sint,cost,phi,densd,idd)
            id=dustarr(ir,it,1)
            taud=taud+kapd(id)*(densd)*dr
            taue=taue+kapd(id)*(dense)*dr
         end do
         write(15,*) sngl(thet*rad2deg),sngl(taue)
     1        ,sngl(taud),sngl(taue+taud)
      end do
      close(15)

c     integrate Av along equatorial direction
      taud=0.d0
      taue=0.d0
      thet=90.d0*deg2rad
      cost=0.d0
      sint=1.d0
      phi=0.d0
      do ir=1,nrg-1
         dr=rarr(ir+1)-rarr(ir)
         rad=0.5d0*(rarr(ir)+rarr(ir+1))
c         call densenv(rad,thet,cost,sint,phi,dense,radamb)
         call densdisk(rad,sint,cost,phi,densd,idd)
         taud=taud+kappav(idd)*rmincgs*(densd)*dr
c     taue=taue+kappav*rmincgs*(dense)*dr
         tauarr(ir)=taud
      end do

c     calculate Av along ethet,ephi directions (peeled image)
      open(unit=15,file='A_v.dat',status='unknown')
      write(15,*) 
     1 'theta   Av_env       Av_disk      Av_tot'
      thet=thete
      cost=cos(thet)
      sint=sin(thet)
      call locate(thetarr,ntg,thet,it)
      phi=0.d0
      taud=0.d0
      taue=0.d0
      do ir=1,nrg-1
         dr=rarr(ir+1)-rarr(ir)
         rad=0.5d0*(rarr(ir)+rarr(ir+1))
         id=dustarr(ir,it,1)
         call densenv(rad,thet,cost,sint,phi,dense,
     1	radamb,ide)
         call densdisk(rad,sint,cost,phi,densd,idd)
         taud=taud+kappav(id)*rmincgs*(densd)*dr
         taue=taue+kappav(id)*rmincgs*(dense)*dr
      end do
      write(15,*) sngl(thet*rad2deg),sngl(taue*1.086)
     1     ,sngl(taud*1.086),sngl((taue+taud)*1.086)
      close(15)

      r=1.d0*autors
      call locate(rarr,nrg,dble(r),ir)
      dens1=densarr(ir,1,1)
      print*,'r,dens1',r,dens1

c     write out density array to read in IDL
c     make the cartesian grid here (fortran is faster than IDL)
c     this makes a grid in the x-z plane
      n=200
c     xmax=sngl(rarr(nrg))
c     xmax=rmaxd
      xmax=sngl(xymaxdens)
      print*,'xmax',xmax
      dx=xmax/real(n)
      do ix=1,n
         x=dx*(real(ix)-0.5)
      do iy=1,n
         y=dx*(real(iy)-0.5)
         if (iy.eq.1.and.ix.eq.1) print*,'x0,y0',x,y
         r=sqrt(x**2+y**2)
         thet=atan2(y,x)
         if (r.gt.rarr(1).and.r.lt.rarr(nrg)) then
            call locate(rarr,nrg,dble(r),ir)
            if (ntg.gt.1) then
               call locate(thetarr,ntg,(thet),it)
            else
               it=1
            end if
            xyarr(ix,iy)=sngl(densarr(ir,it,1))
         else
            xyarr(ix,iy)=0.
         end if
      end do
      end do
      open(unit=15,file='densxz.dat',status='unknown')
      do ix=1,n
         write(15,*) (xyarr(ix,iy),iy=1,n)
      end do
      close(15)

c     write out array in x-y plane at z=.5*rmax to see phi dependence
      thet=pi/2.d0
      if (ntg.gt.1) then
         call locate(thetarr,ntg,(thet),it)
c     print*,'it',it
      else
         it=1
      end if
      do ix=1,n
         x=dx*(real(ix)-0.5)
      do iy=1,n
         y=dx*(real(iy)-0.5)
         if (iy.eq.1.and.ix.eq.1) print*,'x0,y0',x,y
         r=sqrt(x**2+y**2)
         phi=atan2(y,x)
         call locate(rarr,nrg,dble(r),ir)
c     print*,'ir',ir
         if (npg.gt.1) then
            call locate(phiarr,npg,(phi),ip)
         else
            ip=1
         end if
         xyarr(ix,iy)=sngl(densarr(ir,it,ip))
      end do
      end do
      open(unit=15,file='densxy.dat',status='unknown')
      do ix=1,n
         write(15,*) (xyarr(ix,iy),iy=1,n)
      end do
      close(15)

c     write out massarr
      open(unit=15,file='marr.unf',status='unknown',form='unformatted')
      write(15) nrg,ntg,npg
      write(15) (rarr(ir)/autors,ir=1,nrg)
      write(15) (thetarr(it),it=1,ntg)
      write(15) (phiarr(ip),ip=1,npg)
      write(15) (((massarr(ir,it,ip),ir=1,nrg),it=1,ntg)
     1     ,ip=1,npg)
      close(15)

c     write out densarr
      open(unit=15,file='darr.unf',status='unknown',form='unformatted')
      write(15) nrg,ntg,npg
      write(15) (rarr(ir)/autors,ir=1,nrg)
      write(15) (thetarr(it),it=1,ntg)
      write(15) (phiarr(ip),ip=1,npg)
      write(15) (((densarr(ir,it,ip),ir=1,nrg),it=1,ntg)
     1     ,ip=1,npg)
      close(15)
      
c	write out dustarr
      open(unit=15,file='dust.unf',status='unknown',form='unformatted')
      write(15) nrg,ntg,npg
      write(15) (rarr(ir)/autors,ir=1,nrg)
      write(15) (thetarr(it),it=1,ntg)
      write(15) (phiarr(ip),ip=1,npg)
      write(15) (((dustarr(ir,it,ip),ir=1,nrg),it=1,ntg)
     1     ,ip=1,npg)
      close(15)

      end


      







