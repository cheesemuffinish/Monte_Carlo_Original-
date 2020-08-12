      subroutine envset

c     2001/02/17 (baw) constants for envelope assuming TSC
c     sets cavity shape parameters also.

      use tts_mod
      use opacin_mod
      implicit none


      real*8 pi,const,sphmass,rd,rdcgs,rmincgs,z1max,r1max,rtmp,r2max

      pi=4.d0*datan(1.d0)

c     infall;  Rd is disk radius, units of rmin 
c     BAW 7/5/98 make rd = rchole = rc
c     rd = rchole
c     BAW 12/15/99 make rd=rc
      rd=rc
      rdcgs=rd*rstar*rsol
c     rhoe0 is factor in front of density distribution, given by
c     infall calculation.
c     rhoe0=(infallrate)/(4*pi)/sqrt(G*Mcore)/rd**1.5
      const=1.d0/3.15576d7*sqrt(msol)/(4.d0*pi)/sqrt(gn)
      write(6,*)'const',const
      rhoe0=const*rate/sqrt(massc)/rdcgs**1.5
      sphmass=2.d0/3.d0/3.15576d7/sqrt(msol)*rate/sqrt(massc)/
     1     sqrt(2.d0*gn)
      sphmass=sphmass*(rmax*rstar*rsol)**1.5
      write(6,*) 'rhoe0 of envelope',rhoe0
      rmincgs=rstar*rsol
      write(6,*) 'rmincgs',rmincgs
      windmu0=cos(thetmu0*deg2rad)

c     ambient density
      rhoamb=5.0d-20

c     hole in bubble
c     roa=rmax*tand(buboa)
      cosbuboa=cos(buboa*deg2rad)
      
c     1999
c     outflow boundaries
c     z1max is height where opening angle is measured.  take as outer
c     bound.
      z1max=rmax
      r1max=z1max*tan(thet1*deg2rad)
c     z=a+b*x**beta
c     c1=b
c     z01=a
c     ex1=beta
c     r=x  (r is cylindrical radius)
      c1e=(z1max-z01)/r1max**ex1
      print*,'z1max,z01,r1max,ex1,c1e'
      print*,z1max,z01,r1max,ex1,c1e
c     z01 is input
c     check
      if(z01.lt.0.d0) then
         rtmp=(-z01/c1e)**(1.d0/ex1)
         print*,'hole 1 intersects disk at ',rtmp/autors, ' AU'
      end if
      if(ipoly.eq.1) then
         r2max=z1max*tan(thet2*deg2rad)
         c2e=(z1max-z01)/r2max**ex2
         print*,'r2max,ex2,c2e',r2max,ex2,c2e
c     z02 is input
         if(z02.lt.0.d0) then
            rtmp=(-z02/c2e)**(1.d0/ex2)
            print*,'hole 2 intersects disk at ',rtmp/autors, ' AU'
         end if
      end if
      
      return
      end
