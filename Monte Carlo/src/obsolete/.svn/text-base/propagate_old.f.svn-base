c     *********************************************************

      subroutine propagate(iphot,sini,xmax,dum,nphot,idust2
     $     ,icount,iter)

      use tts_mod
      use grid_mod
      use stokesvar_mod
      use taunum_mod
      use opacin_mod
      use vger_mod
      use tauint_mod
      use peeloff_mod
      use dust_mod
      use random
      implicit none

      real*8 xran,xpold,ypold,rpold,rho2,zpold,x,y,xmax
     1   ,sipold,tsum,pol2,sipnew,mu,sth,sux,suy,suz,t4d,t4s
     1   ,dA,t4,t,t4eff,cosbnew,sinbnew,rnew
      real*8 sini(nmu)
      real ran2,dusttemp,dbbfreq,chiR,tnew,told,tave,bbfreq,lucydustfreq
      integer ii(3),dum(3),iphot,iscat,k,ia,ix,iy
     1     ,icount,i,ir,it,ip,inub,iint,iintmax,nphot,ntgh,nacc
     1     ,iabs,idust2,id,iter
     
      integer :: indx(3),i_inter
      
      integer peelid

      iintmax=10000
c      iintmax=100000000
      xpold=xp
      ypold=yp
c     rpold=rp
      zpold=zp
c     already calculated rsq and rtot before call to this routine.

      do i=1,3
         ii(i)=dum(i)
c     print*,'ii ',ii(i)
      end do
            

c     iflag is set to 1 if photon scatters
      iflag=0
      exitflag=0
      aflag=0
      iscat=0
      iabs=0
      iint=0

c     sample tau
      xran=ran()
      tau=-log(xran)	
c     integrate over distance until the optical depth equals tau.
c      print*,'hi,rsq.rtot,zp',rsq,rtot,zp,ii(1),ii(2),ii(3)
      call tauint(iphot,tsum,ii(1),ii(2),ii(3))
c      print*,'hitau,rsq.rtot,zp',rsq,rtot,zp,ii(1),ii(2),ii(3)
      ir=ii(1)
      it=ii(2)
      ip=ii(3)
      idust2=dustarr(ir,it,ip)
      if(exitflag.eq.1) go to 300
      if(aflag.eq.1) then
         write(6,*) 'shouldnt be here, aflag=1'
         go to 400
      end if
c      print*,'after tauint,rtot,ir',rtot,ir
c      call testgrid(xp,yp,zp,rsq,rtot,ii,
c     +     'after tauint',iphot,r2p,cost)
      tot=tot+1.d0
      
c     photon scatters until exit exits disk
c      do while(iint.lt.iintmax)
 30   continue
c     sip=sip*rlam
c     sqp=sqp*rlam
c     sup=sup*rlam
c     svp=svp*rlam
      xran=ran()
c     print*,'rlam',rlam
      if(xran.le.albedo(idust2)) then
         iflag=1
         iscat=iscat+1
c         if (iscat.gt.10000) then
c            print*,'iscat big',iscat
c            print*,'albedo,wave',albedo,wave
c         end if
         iint=iint+1
c     if (sip.lt.(1.0d-3*rlam)) then
c     aflux=aflux+1
c     go to 400
c     end if
         sipnew=sip
         if (ipeel.eq.1) then
           do peelid=1,npeel
             call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp
     1           ,cost,sint,cosp,sinp,phi
     1           ,hit,htot,rsq,rtot
     1           ,tsum,ii(1),ii(2),ii(3),idust,iflag,iphot,peelid)
           end do
         end if
c     NOTE**** need to write vger_3d.f subroutine*******
c         if(iveeg.eq.1) then
c            call vger(xp,yp,zp,sip,sqp,sup,svp
c     1           ,cost,sint,cosp,sinp,phi
c     1           ,pi,r2p,hit,htot,rsq,rtot
c     1           ,tsum,ii,idust,iflag,iphot,eps)
c         end if
         
         pol2=sqp**2+sup**2+svp**2
         if (pol2 .gt. sip**2) then
            print*,'error, P^2, I^2 ', pol2,sip**2
            continue
         end if      
   
c     xran=ran()
c     if(xran.le.rlam) then
         if (idust.lt.2) then
            call stokes(idust2)
         else
c             call stokes6(idust2)
		print*,'oops, error, idust is wrong'
		stop
         end if
         pol2=sqp**2+sup**2+svp**2
         if (pol2 .gt. sip**2) then
            print*,'error, P^2, I^2 ', pol2,sip**2
            continue
         end if
c     iscat=iscat+1
      else
c     aflux=aflux+1
c     go to 400
c     photon absorbed, reemitted at longer wavelengths 
c     set photon origin to disk or envelope
         if (idust2.eq.1.or.idust2.eq.2) then
            i_orig=2
         else
            i_orig=3
         end if
c         i_orig=2
c     iabs=iabs+1
c         if (iabs.gt.1000) then
c            print*,'iabs big',iabs
c         end if
         if(albedo(idust2).eq.1.d0)print*,
     1	'error! albedo=1, yet photon absorbed'
c     print*,'freq before abs',nub
         if (.not.diffus(ir,it,ip)) then
            nabs(ir,it,ip)=nabs(ir,it,ip)+1
            if (ilucy.eq.0) then
               told=tdust(ir,it,ip)
               tdust(ir,it,ip)=dusttemp(ir,it,ip,dble(nabs(ir,it,ip))
     1           ,sngl(tdust(ir,it,ip)),rstar,tstarave,nphot,idust2)
               tnew=tdust(ir,it,ip)
               if (nabs(ir,it,ip).lt.2) then
                  tave=0.25*told+0.75*tnew
               else
                  tave=tnew
               end if
            else
               tave=tdust(ir,it,ip)
c               print*,tave*11605.d0
            end if
            call reemit(tave,idust2,iter)
            uy=sint*sinp
            ux=sint*cosp
            uz=cost
c            if (ir.eq.10.and.it.eq.5.and.ip.eq.1) then
c               print*,'told,tnew,nabs,tave',told*11605.d0,tnew*11605.d0
c     $              ,nabs(ir,it,ip),tave*11605.d0,iter
cc     $              ,nabs(ir,it,ip),nphot,rstar,tstar*11605.d0,idust2
c            end if

         else
c     nabs in diffusion layer only consists of accretion photons
c     march through theta array to find theta of first non-diffusion
c     cell
c     call emitdiff(ir,it,ip,xp,yp,zp)
c     print*,'in diffusion region'

c     find diffusion direction
            if (abs(diffdir(ir,it,ip)) .eq. 1) then
c            print*,'r diffusion'
               if (diffdir(ir,it,ip).eq. 1) then
                 write(*,*) 'error, positive r diffusion not supported'
                  stop
               end if

c              diffuse to front of disk

c              print*,'check: rtot,rarrs',rtot,rarr(ir),rarr(ir+1)
c              print*,'rsq,arrays',rsq,r2arr(ir),r2arr(ir+1)
c              print*,'diffdir',diffdir(ir,it,ip)
c               print*,'should be same',rtot,ir
c               print*,ii(1),ir
c               call testgrid(xp,yp,zp,rsq,rtot,ii,
c     +              'before diffusion',iphot,r2p,cost)
c     okay, we are where we think we are

               do while (diffdir(ir,it,ip) .eq. -1) !find diff surface
                  ir=ir-1
               end do
               ir=ir+1 !point to lower wall of 1st diffusion cell
               rnew=rarr(ir)
               on_wall = .true.
c               if (rnew.gt.rtot) print*,'rnew>rtot',rnew,rtot
c               print*,'eps,rnew,rtot',eps,rnew,rtot
c               rnew=0.9999999d0*rarr(ir) !move photon into cell above diff layer
               xp=xp*rnew/rtot
               yp=yp*rnew/rtot
               rp=sqrt(xp**2+yp**2)
               zp=zp*rnew/rtot
c               rsq=zp**2+xp**2+yp**2
c               rtot=sqrt(rsq)
               rtot=rnew
               rsq=rnew*rnew
               
               sux=-xp/rtot  !unit vector perp to diffusion surf
               suy=-yp/rtot
               suz=-zp/rtot

               call isotrp(sux,suy,suz,ux,uy,uz) !emit isotropically from surf
               cost=uz
               sint=sqrt(1.d0-cost*cost)
               cosp=ux/sint
               sinp=uy/sint
               phi=atan2(sinp,cosp)

c              find new disk surface temperature

               dA=2.d0*pi*rarr(ir)**2*(costarr(it)-costarr(it+1)) !assumes 2-D

               T4d=tdust(ir,it,ip)*tdust(ir,it,ip)
               T4d=T4d*T4d
               
               T4s=Tstar*Tstar
               T4s=T4s*T4s
               
               T4=T4d+4*pi*T4s/(nphot*dA) !surf temp
               Tdust(ir,it,ip)=min(1600.d0/11605.d0,T4**0.25)

c              determine temperature in disk interior

               nacc=0
               i=ir
               do while (diffdir(i,it,ip) .eq. -1)
                  nacc=nacc+nabs(i,it,ip) !accumulate accretion photon counter
                  i=i+1
               end do

               T4d=T4 !initialize integration variable to surf temp
               i=ir+1
               do while (diffdir(i,it,ip) .eq. -1)
                     idust2=dustarr(i,it,ip)
                     T4d=T4d+(3.*pi*T4s/(nphot*dA))*nacc
     +                 *kappav(idust2)*rsol*rstar*densarr(i-1,it,ip)
     +                *chiR(sngl(tdust(i-1,it,ip)),idust2)*(rarr(i)-
     +		    rarr(i-1))
                     tdust(i,it,ip)=min(1600.d0/11605.d0,T4d**0.25)
                     nacc=nacc-nabs(i-1,it,ip)
                     i=i+1
               end do

c              get photon temperature using grey atm approximation

               T4eff=4.*pi*T4s*nacc/(dA*nphot)
               
               mu=ux*sux+uy*suy+uz*suz
               T=min((T4+0.75*mu*T4eff)**0.25,1600.d0/11605.d0)
               
c               nub=dBBfreq(T,idust2)
               if (ilucy.eq.0) then
                  nub=dBBfreq(T,idust2)
               else
c                  nub=BBfreq(T)
c                 nub=lucydustfreq(T,idust2)
                  nub=0.
               end if
               
               ir=ir-1 !move pointer to first cell below diffusion layer
               
               ii(1)=ir
               
               
            else if (abs(diffdir(ir,it,ip)) .eq. 2) then

c              diffuse to top or bottom of disk

               if (diffdir(ir,it,ip) .eq. -2) then
                  do while (diffus(ir,it,ip))
                     it=it-1
                  end do
                  it=it+1
               else
                  do while (diffus(ir,it,ip))
                     it=it+1
                  end do
               end if
c     print*,'r,it',rarr(ir),it
c     get photon position at diffusion layer
               cosb=zp/rtot
               sinb=sqrt(1.d0-cosb*cosb)
               cosbnew=1.00001*costarr(it)
               sinbnew=sqrt(1.d0-cosbnew*cosbnew)
               xp=xp*sinbnew/sinb !move photon
               yp=yp*sinbnew/sinb
               rp=sqrt(xp**2+yp**2)
               zp=rtot*cosbnew
c     zp=rtot*sign(1.00001*abs(costarr(it)),costarr(it))
c     calculate mu of disk at position of photon
               mu=costarr(it)
               sth=sintarr(it)
c     unit vector perpendicular to surface
               if (mu .gt. 0.d0) then
                  sux=-mu*xp/sth/rtot
                  suy=-mu*yp/sth/rtot
                  suz=sth
               else
                  sux=mu*xp/sth/rtot
                  suy=mu*yp/sth/rtot
                  suz=-sth
               end if
c     get isotropic emission for photon
               call isotrp(sux,suy,suz,ux,uy,uz)
               cost=uz
               sint=sqrt(1.d0-cost*cost)
               cosp=ux/sint
               sinp=uy/sint
               phi=atan2(sinp,cosp)
c     reassign cost,cosp, etc
               if (costarr(it).lt.0.d0) then
                  it=it-1
               end if

c              find new disk surface temperature

               dA=pi*sth*(rarr(ir+1)**2-rarr(ir)**2) !assumes 2-D

               T4d=tdust(ir,it,ip)*tdust(ir,it,ip)
               T4d=T4d*T4d
               
               T4s=Tstar*Tstar
               T4s=T4s*T4s
               
c     T4=T4d+4*pi*Rstar**2*T4s/(Lsfrac*Nphot*dA)
               T4=T4d+4*pi*T4s/(nphot*dA) !surf temp
               Tdust(ir,it,ip)=min(1600.d0/11605.d0,T4**0.25)

c              determine temperature in disk interior

               T4d=T4 !initialize integration variable to surf temp

               if (diffdir(ir,it,ip) .eq. -2) then

                  nacc=0
                  i=it
                  do while (diffdir(ir,i,ip) .eq. -2)
                     nacc=nacc+nabs(ir,i,ip) !accumulate accretion photon counter
                     i=i+1
                  end do

                  i=it+1
                  do while (diffdir(ir,i,ip) .eq. -2)
                     idust2=dustarr(ir,i,ip)
                     T4d=T4d+(3.*pi*T4s/(nphot*dA))*nacc
     +                 *kappav(idust2)*rsol*rstar*densarr(ir,i-1,ip)
     +                    *chiR(sngl(Tdust(ir,i-1,ip)),idust2)
     +                    *0.5d0*(rarr(ir)+rarr(ir+1))
     +                    *(thetarr(i)-thetarr(i-1))
                     tdust(ir,i,ip)=min(1600.d0/11605.d0,T4d**0.25)
                     nacc=nacc-nabs(ir,i-1,ip)
                     i=i+1
                  end do

               else

                  nacc=0
                  i=it
                  do while (diffdir(ir,i,ip) .eq. 2)
                     nacc=nacc+nabs(ir,i,ip) !accumulate accretion photon counter
                     i=i-1
                  end do

                  i=it-1
                  do while (diffdir(ir,i,ip) .eq. 2)
                     idust2=dustarr(ir,i+1,ip)
c                     idust2=dustarr(ir,i,ip)
                     T4d=T4d+(3.*pi*T4s/(nphot*dA))*nacc
     +                 *kappav(idust2)*rsol*rstar*densarr(ir,i+1,ip)
     +                    *chiR(sngl(tdust(ir,i+1,ip)),idust2)
     +                    *0.5d0*(rarr(ir)+rarr(ir+1))
     +                    *(thetarr(i+2)-thetarr(i+1))
                     tdust(ir,i,ip)=min(1600.d0/11605.d0,T4d**0.25)
                     nacc=nacc-nabs(ir,i+1,ip)
                     i=i-1
                  end do


               end if

c              get photon temperature using grey atm approximation

c     T4eff=3.*Lacc*(Rstar/r)**3*(1.d0-sqrt(Rstar/r))*T4s
               T4eff=4.*pi*T4s*nacc/(dA*nphot)
               
               mu=ux*sux+uy*suy+uz*suz
               T=min((T4+0.75*mu*T4eff)**0.25,1600.d0/11605.d0)

c               nub=dBBfreq(T,idust2)
               if (ilucy.eq.0) then
                  nub=dBBfreq(T,idust2)
              else
c                  nub=BBfreq(T)
c                 nub=lucydustfreq(T,idust2)
                 nub=0.
               end if               

               if (diffdir(ir,it,ip) .eq. -2) then
                  it=it-1
               else
                  it=it+1
               end if
               
               ii(2)=it
               
            else if (abs(diffdir(ir,it,ip)).eq.3) then
               write(*,*) 'error, phi diffusion not supported'
               stop
            else
               write(*,*) 'error, unknown diffusion direction'
               stop
            end if

         end if

c         call testgrid(xp,yp,zp,rsq,rtot,ii,
c     +              'after diffusion',iphot,r2p,cost)

c         nabs(ir,it,ip)=nabs(ir,it,ip)+1
c         tdust(ir,it,ip)=dusttemp(ir,it,ip,nabs(ir,it,ip)
c     1        ,tdust(ir,it,ip),rstar,tstarave,nphot)
c         call reemit(sngl(tdust(ir,it,ip)))

         iflag=0
c     print*,'freq after abs',nub
         call opacset(nub)
         wave=1.2398d0/nub
c         if (wave.lt.0.8) then 
c            print*,'thermal photon with wave = ',wave,'um'   
c            print*,'sip',sip
c         end if
c     wave=2.99792458d14/nub
	   idust2=dustarr(ir,it,ip)
	   do id=1,ndg
            kapd(id)=kappa(id)*rstar*rsol*kappav(id)
         end do
         sipnew=sip
         if (ipeel.eq.1) then
           do peelid=1,npeel
            call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp
     1           ,cost,sint,cosp,sinp,phi
     1           ,hit,htot,rsq,rtot
     1           ,tsum,ii(1),ii(2),ii(3),idust,iflag,iphot,peelid)
           end do
         end if
c         if(iveeg.eq.1) then
c            call vger(xp,yp,zp,sip,sqp,sup,svp
c     1           ,cost,sint,cosp,sinp,phi
c     1           ,pi,r2p,hit,htot,rsq,rtot
c     1           ,tsum,ii,idust,iflag,iphot,eps)
c         end if
         
c     aflux=aflux+1.d0
         iint=iint+1
      end if
c     sample tau
      xran=ran()
      tau=-log(xran)
c     integrate distance until optical depth equals tau.
      xpold=xp
      ypold=yp
      zpold=zp
c     not sure if we need to do this here, might already have this
c     rsq=xp**2+yp**2+zp**2
c     rtot=sqrt(rsq)
c     rpold=rp
      ux=sint*cosp
      uy=sint*sinp
      uz=cost

c      print*,'end,rsq.rtot,zp',rsq,rtot,zp,ii(1),ii(2),ii(3)
      call tauint(iphot,tsum,ii(1),ii(2),ii(3))
      ir=ii(1)
      it=ii(2)
      ip=ii(3)
      idust2=dustarr(ir,it,ip)
      if(exitflag.eq.1) go to 300
c      print*,'after tauint,rtot,ir',rtot,ir
c      call testgrid(xp,yp,zp,rsq,rtot,ii,
c     +     'after tauint',iphot,r2p,cost)

c     no longer absorbing photon by star, reemitting
c      if (aflag.eq.1) then
c         write(6,*) 'photon absorbed by star, iphot',iphot
c         aflux=aflux+1
c     reemit thermal photon
c         go to 400
c      end if 

c     testing!!!!!!
c      if (iint.gt.iintmax) then
c         abflux=abflux+1
c         go to 400
c      end if

      go to 30

c     should never be here
c     print*,' oops, should not be here, propagate, L145'

 300  continue
c     photon exits
c      if (iint.gt.10000) then
c         print*,'large number of interactions',iint
c      end if

c     bin angle
c     NOTE:  assuming axisymmetric, and z=-z.  
c      k=int(real(nmu-1)*abs(cost)+1.5d0)
c     for comparison to jon
c      k=min(int(real(nmu)*abs(cost)+1),nmu)
c     20080826 new binning in cost, not combining about midplane
c      print*,k
      k=int(real(nmu)*(1-cost)/2.)+1
      if(k.lt.0.or.k.gt.nmu) then
         print*,'cost binning error, k, nmu',k,nmu
         print*, 'cost,sint',cost,sint
         print*, 'iphot',iphot
      end if      

      if (phi.lt.0.d0) phi=phi+r2p
      ip=int(nph*phi/r2p)+1
      if(ip.lt.0.or.ip.gt.nph) then
         print*, 'phi binning error, m, nph ',ip,nph
         print*, 'phi ', phi
         print*, 'iphot ',iphot
      end if

c     old imaging.  (not even used anymore but don't ever want to have
c     to figure these out again)
c      if (cost.lt.0.d0) then
c         x = -ypold*cosp+xpold*sinp
c         y = -zpold*sini(k)-ypold*u(k)*sinp
c     &        -xpold*u(k)*cosp
c      else
c         x = ypold*cosp-xpold*sinp
c         y = zpold*sini(k)-ypold*u(k)*sinp
c     &        -xpold*u(k)*cosp
c      end if
 
c     20080826 imaging with no mirror of cost (not used anymore)
      x = ypold*cosp-xpold*sinp
      y = zpold*sini(k)-ypold*u(k)*sinp
     &     -xpold*u(k)*cosp
     
c     first, sum fluxes
      rho2=(x**2+y**2)
      if (rho2.lt.aperture2(nap)*1.0001d0) then
         flux=flux+1.
         if (iflag.eq.1) then
            sflux=sflux+1
            nscat=nscat+iscat
         end if
c      else
c         print*,'who am I?',iphot
c     print*,'rho2 bigger than aperture2(nap)',rho2,aperture2(nap)
c     well, maybe the user wants rho2 bigger than aperture2(nap)
      end if

c     find frequency bin, and make sure it is inside the limits
      inub=min(max(int(nfreq*log(nub/numin)/lnurat)+1,1),nfreq)
      if (inub.lt.1.or.inub.gt.nfreq) then
         print*,'inub out of limits!,inub,nub',inub,nub
      end if

      if (i_orig.lt.1.or.i_orig.gt.3) print*,'error, IOR wrong',i_orig

c     loop over apertures, and bin photon into SEDs

      if(iflag==1) then
        i_inter=2
      else
        if(i_orig==1) then
          i_inter=1
        else
          i_inter=3
        end if
      end if

      indx = (/1,1+i_orig,5+i_inter/)
      !    this is an array with 3 indices,  the first is 1, the second is the origin of the photon, the third is the interaction type
      !     1 = all
      !     2 = star origin
      !     3 = disk origin (either accretion or reeimision)
      !     4 = envelope origin (reemission)
      !     5 = external illumination origin
      !     6 = direct star 
      !     7 = scattered
      !     8 = thermal emission
      

      do ia=1,nap
        if(rho2 < aperture2(ia)) then

c         add photon to flux arrays
          si(inub,k,ip,ia,indx)=si(inub,k,ip,ia,indx)+sip
          sq(inub,k,ip,ia,indx)=sq(inub,k,ip,ia,indx)+sqp
          su(inub,k,ip,ia,indx)=su(inub,k,ip,ia,indx)+sup
          sv(inub,k,ip,ia,indx)=sv(inub,k,ip,ia,indx)+svp
      
c         add photon to flux error arrays
          si2(inub,k,ip,ia,indx)=si2(inub,k,ip,ia,indx)+sip*sip
          sq2(inub,k,ip,ia,indx)=sq2(inub,k,ip,ia,indx)+sqp*sqp
          su2(inub,k,ip,ia,indx)=su2(inub,k,ip,ia,indx)+sup*sup
          sv2(inub,k,ip,ia,indx)=sv2(inub,k,ip,ia,indx)+svp*svp
      
          nums(inub,k,ip,ia,indx)=nums(inub,k,ip,ia,indx)+1.d0
      
          aveinc(inub,k,ip,ia,indx)=aveinc(inub,k,ip,ia,indx)+abs(cost)
    
        end if
      end do
      
c      if (iflag.eq.1.and.i_orig.eq.1) scount=scount+1.d0
      if (iflag.eq.1) scount=scount+1.d0
      
c     testing**************
c    removed images except for peeling off
c     if (cost.ge.0) then      
c      ix=int((x+xmax)/fractx)+1
c      iy=int((y+xmax)/fractx)+1
c      if (ix.le.nx.and.iy.le.nx.and.ix.gt.0.d0.and.iy.gt.0.d0) then
c         image(ix,iy,k)=image(ix,iy,k)+(sip)
c         imagei(ix,iy,k)=imagei(ix,iy,k)+(sip)
c         imageq(ix,iy,k)=imageq(ix,iy,k)+(sqp)
c         imageu(ix,iy,k)=imageu(ix,iy,k)+(sup)
c         imagev(ix,iy,k)=imagev(ix,iy,k)+(svp)
c         image2(ix,iy,k)=image2(ix,iy,k)+sip*sip
c         numi(ix,iy,k)=numi(ix,iy,k)+1
c      else
c         icount=icount+1
c      end if
      
 400  continue

      return
      end
      
      
c     *********************************************************



