c     *********************************************************

      subroutine output(iwav,cwav,ctherm,MXWAV)
      
      use tts_mod
      use grid_mod
      use stokesvar_mod
      use opacin_mod
      use vger_mod
      use out_mod
      implicit none 


      integer iwav,MXWAV
      character cwav(MXWAV)*(*),ctherm(MXWAV)*(*)
      character beg*9

      real*8 tmpsi,dflux,scatave,cpusec,cpuhrs,cpumin,sec
     1     ,fq,etime,fqp,fup,fvp,fi,fu,fv,rnorm,fp,totvg,fip
     1     ,rad,thet,x,z,tave1,tave2,count

      real time(2)
      
      real*8 fig(napmax)
      character fchar*3

      integer i,ia,imin,ihrs,ix,iy,it,ip,inub,ir,ift
      
      integer MXNAM
      parameter (MXNAM=20)
      character filnam(MXNAM)*25
      
c     print*,'iwave,cwave',iwav,cwav(iwav)
      call namer(1,cwav(iwav),filnam,MXNAM)

      open(unit=21,file=filnam(8),status='unknown')
      do it=1,nmu
      do inub=1,nfreq
         wave=1.2398d0/(numin*nurat**((real(inub)-0.5)/nfreq))
         write(21,*)wave, (star(inub,it)/nunorm)
      end do
      end do
      close(21)

c     flux.dat and sflux.dat
      open(unit=11,file=filnam(6),status='unknown')
      open(unit=13,file=filnam(7),status='unknown')
      write(11,*) 
     1'wave,          i,          q/i,          u/i,          v/i'
      write(13,*) 
     1'wave,          i,          q/i,          u/i,          v/i'
      rnorm=real(np)/real(nmu)
      do ia=1,nap
      do i=1,nmu
      do inub=1,nfreq
         wave=1.2398d0/(numin*nurat**((real(inub)-0.5)/nfreq))
         tmpsi=si(inub,i,ia,1)+si(inub,i,ia,2)+si(inub,i,ia,3)
         if(spi(inub,i,ia).eq.0.0d+00) then
            fip=0.
            fqp=0.
            fup=0.
            fvp=0.
         else
            fip=(spi(inub,i,ia)/nunorm)
            fqp=(spq(inub,i,ia)/spi(inub,i,ia))
            fup=(spu(inub,i,ia)/spi(inub,i,ia))
            fvp=(spv(inub,i,ia)/spi(inub,i,ia))
         end if
         if(tmpsi.eq.0.d0) then
            fi=0.
            fq=0.
            fu=0.
            fv=0.
         else
            fi=(tmpsi/nunorm)
            fq=(sq(inub,i,ia)/tmpsi)
            fu=(su(inub,i,ia)/tmpsi)
            fv=(sv(inub,i,ia)/tmpsi)
         end if
         fip=(spi(inub,i,ia)/rnorm/nunorm)
         write(11,*) (wave),fi,fq,fu,fv
c     1 ,(u(i))
         write(13,*) (wave),fip,fqp,fup,fvp
c     1 ,(u(i))
      end do
      end do
      end do
      close(11)
      close(13)

c     Write flux in form suitable for grid - Added 28th June 2005

      write(fchar,'(I3)') nap
      open(unit=11,file='flux_grid.dat',status='unknown')
      write(11,'(2X,"LAMBDA",2X,'//fchar//'(4X,F8.0))') 
     &(sqrt(aperture2(ia))/autors,ia=1,nap)
      do i=1,nmu
         do inub=1,nfreq
         wave=1.2398d0/(numin*nurat**((real(inub)-0.5)/nfreq))
         do ia=1,nap
            tmpsi=si(inub,i,ia,1)+si(inub,i,ia,2)+si(inub,i,ia,3)
            if(tmpsi.eq.0.d0) then
               fig(ia)=0.
            else
               fig(ia)=(tmpsi/nunorm)
            end if
         end do
         write(11,'(F10.4,'//fchar//'(2X,ES10.3))')
     1 (wave),(fig(ia),ia=1,nap)
      end do
      end do
      close(11)



c     fluxdisk.dat, fluxstar.dat, fluxenv.dat, origin of photons.
c     didn't track all four stokes vectors for these so write out zeros for
c     q,u,v
      open(unit=11,file='fluxstar.dat',status='unknown')
      open(unit=13,file='fluxdisk.dat',status='unknown')
      open(unit=14,file='fluxenv.dat',status='unknown')
      write(11,*) 
     1'wave,          i,          q/i,          u/i,          v/i'
      write(13,*) 
     1'wave,          i,          q/i,          u/i,          v/i'
      write(14,*) 
     1'wave,          i,          q/i,          u/i,          v/i'
      rnorm=real(np)/real(nmu)
      do ia=1,nap
      do i=1,nmu
      do inub=1,nfreq
         wave=1.2398d0/(numin*nurat**((real(inub)-0.5)/nfreq))
         tmpsi=si(inub,i,ia,1)
         if(tmpsi.eq.0.d0) then
            fi=0.
         else
            fi=(tmpsi/nunorm)
         end if
         fq=0.
         fu=0.
         fv=0.
         write(11,*) (wave),fi,fq,fu,fv
         tmpsi=si(inub,i,ia,2)
         if(tmpsi.eq.0.d0) then
            fi=0.
         else
            fi=(tmpsi/nunorm)
         end if
         write(13,*) (wave),fi,fq,fu,fv
         tmpsi=si(inub,i,ia,3)
         if(tmpsi.eq.0.d0) then
            fi=0.
         else
            fi=(tmpsi/nunorm)
         end if
         write(14,*) (wave),fi,fq,fu,fv
      end do
      end do
      end do


      close(11)
      close(13)
      close(14)

      open(unit=11,file=filnam(9),status='unknown')
	write(11,*) 'standard deviations of i,q,u,v'
	write(11,*) 'normalize to i for q,u,v to compare directly with'
	write(11,*) '    flux file (assumes i has no error) '
      write(11,*) 
     1'  wave,       sig_i,         sig_q/i,      sig_u/i,      sig_v/i'
      do ia=1,nap
      do i=1,nmu
      do inub=1,nfreq
         wave=1.2398d0/(numin*nurat**((real(inub)-0.5)/nfreq))
         tmpsi=si(inub,i,ia,1)+si(inub,i,ia,2)+si(inub,i,ia,2)
         if(nums(inub,i,ia).le.1) then
            fi=0.
            fq=0.
            fu=0.
            fv=0.
            fp=0.
         else
            fi=(si2(inub,i,ia)*tmpsi)/nunorm
            fq=(sq2(inub,i,ia))
            fu=(su2(inub,i,ia))
            fv=(sv2(inub,i,ia))
            fp=1./sqrt(real(nums(inub,i,ia)))*(tmpsi)
	   end if
         write(11,*) (wave),fi,fq,fu,fv
      end do
      end do
      end do
      close(11)


c     flux
c     flux=flux+aflux
c     flux which hits disk
      dflux=(tot/flux)
c     scattered flux
      if(nscat.gt.0) then
         scatave=nscat/sflux
      else
         scatave=0.
      end if
      sflux=sflux/flux
c     absorbed flux
c     aflux=aflux/np
c     abflux=abflux/np
c     cpu time
      cpusec=etime(time)
      cpuhrs=cpusec/3600.
      ihrs=int(cpuhrs)
      cpumin=(cpuhrs-ihrs)*60.
      imin=int(cpumin)
      sec=(cpumin-imin)*60.

      write(12,*) 'rmine used ',rmine
      write(12,*) 'rmax,rmaxd,zmax',rmax,rmaxd,zmax
      write(12,*) 'rmaxi,rmind',rmaxi,rmind
      write(12,*) 'photons   ',np
      write(12,*) 'taur ',taur
      write(12,*) 'massenv, massdisk ',massenv,massdisk
      write(12,*) 'total flux', flux
      write(12,*) 'flux which hits disk+envelope',dflux,dflux*flux
      write(12,*) 'scattered flux',sflux,sflux*flux
      write(12,*) 'ave number of scatters in this flux',scatave
      write(12,*) 'flux which gets absorbed (and reemitted) by star',
     1     aflux,aflux/np
      write(12,*) 'killed photons,flux ',abflux,abflux/np
      write(12,*) 'cputime: ',ihrs,' hrs ',imin,' min ',sec,' sec'

      it=int(thete*180.*1.001/pi+1)
      call namer(it,cwav(iwav),filnam,MXNAM)

      if (ipeel.eq.1) then
         open(unit=13,file='e'//filnam(6),status='unknown')
	   write(13,*) 'flux and polarization in the 3 apertures '
         write(13,*)
     1'  I              Q/I              U/I               V/I'
         write(13,*) 'total flux'
         do i=1,nap
         do inub=1,nfreq
            wave=1.2398d0/(numin*nurat**((real(inub)-0.5)/nfreq))
            if(ti(inub,i).eq.0.d0) then
               fi=0.
               fq=0.
               fu=0.
               fv=0.
            else
               fi=(ti(inub,i)/nunorm)
               fq=(tq(inub,i)/ti(inub,i))
               fu=(tu(inub,i)/ti(inub,i))
               fv=(tv(inub,i)/ti(inub,i))
            end if
            write(13,*) wave,fi,fq,fu,fv
         end do
         end do
         write(13,*) 'scattered flux only'
         do i=1,nap
         do inub=1,nfreq
            if(tis(inub,i).eq.0.d0) then
               fi=0.
               fq=0.
               fu=0.
               fv=0.
            else
               fi=(tis(inub,i)/nunorm)
               fq=(tqs(inub,i)/tis(inub,i))
               fu=(tus(inub,i)/tis(inub,i))
               fv=(tvs(inub,i)/tis(inub,i))
            end if
            write(13,*) wave,fi,fq,fu,fv
         end do
         end do
         write(13,*) 'direct stellar flux'
c         do i=1,nap
         do inub=1,nfreq
            write(13,*) wave,(estar(inub)/nunorm)
         end do
c         end do
	   write(13,*) ' standard deviations in the 3 apertures '
         write(13,*)
     1' sigI          sig(Q/I)     sig(U/I)     sig(V/I)     Poisson(I)'
         do ia=1,nap
         do inub=1,nfreq
            if(numt(inub,ia).le.1) then
               fi=0.
               fq=0.
               fu=0.
               fv=0.
               fp=0.
            else
               fi=(ti2(inub,ia)*ti(inub,ia)/nunorm)
               fq=(tq2(inub,ia))
               fu=(tu2(inub,ia))
               fv=(tv2(inub,ia))
               fp=1./sqrt(real(numt(inub,ia)))*(ti(inub,ia))
            end if
            write(13,*) fi,fq,fu,fv,fp
         end do
         end do
         close(13)

         do ift=1,nbnd
            if (ift.eq.1)  beg='e_1.14um_'
            if (ift.eq.2)  beg='e_1.61um_'
            if (ift.eq.3)  beg='e_2.08um_'
            if (ift.eq.4)  beg='e_3.60um_'
            if (ift.eq.5)  beg='e_4.50um_'
            if (ift.eq.6)  beg='e_5.80um_'
            if (ift.eq.7)  beg='e_8.00um_'
            if (ift.eq.8)  beg='e_24.0um_'
            if (ift.eq.9)  beg='e_70.0um_'
            if (ift.eq.10) beg='e_160.um_'
            if (ift.eq.11) beg='e_7.9.um_'
            if (ift.eq.12) beg='e_8.8.um_'
            if (ift.eq.13) beg='e_9.7.um_'
            if (ift.eq.14) beg='e_10.3.m_'
            if (ift.eq.15) beg='e_11.7.m_'
            if (ift.eq.16) beg='e_12.5.m_'
            if (ift.eq.17) beg='e_17.9.m_'
            if (ift.eq.18) beg='e_20.8.m_'
            if (ift.eq.19) beg='e_24.5.m_'
            if (ift.eq.20) beg='e_443.um_'
            if (ift.eq.21) beg='e_800.um_'
            open(unit=21,file=beg//filnam(2),status='unknown')
            open(unit=22,file=beg//filnam(3),status='unknown')
            open(unit=23,file=beg//filnam(4),status='unknown')
            open(unit=24,file=beg//filnam(5),status='unknown')
            do iy=1,nxhst
               write(21,*) ((image_b_i(ix,iy,ift)),ix=1,nxhst)
               write(22,*) ((image_b_q(ix,iy,ift)),ix=1,nxhst)
               write(23,*) ((image_b_u(ix,iy,ift)),ix=1,nxhst)
               write(24,*) ((image_b_v(ix,iy,ift)),ix=1,nxhst)
            end do
            close(21)
            close(22)
            close(23)
            close(24)
         end do
      end if

      if (iveeg.eq.1) then
         open(unit=21,file=filnam(10),status='unknown')
	   write(21,*) '      F         Q/F         U/F           V/F'
	   totvg=ivgflux
	   write(21,*) (totvg),(qvgflux/totvg),
     1   (uvgflux/totvg),(vvgflux/totvg)
	   write(21,*) '    errF       errQ/F      err(U/F)      err(V/F)'
	   write(21,*) (ivgerr)*(totvg),(qvgerr),
     1   (uvgerr),(vvgerr)
	   close(21)
c        if you really want vger images, have to assign them names
c        in namer.  for now, just assume doing one wavelength.
         open(unit=21,file=filnam(11),status='unknown')
         open(unit=22,file=filnam(12),status='unknown')
         open(unit=23,file=filnam(13),status='unknown')
         open(unit=24,file=filnam(14),status='unknown')
	   do it=1,ncvg
	     write(21,*) (ivg(it,ip),ip=1,npvg)
	     write(22,*) (qvg(it,ip),ip=1,npvg)
	     write(23,*) (uvg(it,ip),ip=1,npvg)
	     write(24,*) (vvg(it,ip),ip=1,npvg)
	   end do
	   close(21)
	   close(22)
	   close(23)
	   close(24)
	end if

c      do ip=1,npg
c      do it=1,ntg
c      do ir=1,nrg
c         if (nabs(ir,it,ip).eq.0) tdust(ir,it,ip)=0.d0
c      end do
c      end do
c      end do

      print*,'rmin',rarr(1)

c     write out tarr.unf in tfinal.f
c     open(unit=15,file='tarr.unf',status='unknown',form='unformatted')
c     write(15) nrg,ntg,npg
c     write(15) (rarr(ir)/autors,ir=1,nrg)
c     write(15) (thetarr(it),it=1,ntg)
c     write(15) (phiarr(ip),ip=1,npg)
c     write(15) (((tdust(ir,it,ip)*11605.d0,ir=1,nrg),it=1,ntg)
c     1     ,ip=1,npg)
c     close(15)

      open(unit=15,file='nabs.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((nabs(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

      open(unit=15,file='tmidplane.dat',status='unknown')
      it=(ntg+1)/2
      write(15,*) 'r(rstar),r(au),tau,Tdust,Tdust2 at thet=',
     $     thetarr(it)*rad2deg
      do ir=1,nrg-1
         write(15,*) ir,(ravearr(ir)),(ravearr(ir)/autors),
     $        (tauarr(ir)),(tdust(ir,it,1))*11605.d0,
     $        (tdust2(ir,it,1))*11605.d0
      end do
      close(15)

      open(unit=15,file='tave.dat',status='unknown')
      write(15,*) 'r(rstar),r(au),  Tave, Tave2'
      do ir=1,nrg-1
         tave1=0.d0
         tave2=0.d0
         count=0.d0
         do it=1,ntg-1
         do ip=1,npg-1
            tave1=tdust(ir,it,ip)+tave1
            tave2=tdust2(ir,it,ip)+tave2
            count=count+1.d0
         end do
         end do
         tave1=tave1/count
         tave2=tave2/count
         write(15,*) ir,(ravearr(ir)),(ravearr(ir)/autors),
     $        tave1*11605.d0,tave2*11605.d0
      end do
      close(15)

      return
      end

