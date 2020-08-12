      subroutine stokespeel(idust,sip,sqp,sup,svp,cost,sint
     1     ,hit,htot,pi,r2p,ri1,bmu,costp,sintp,idust2)
      
c     history:
c     2003/12/26 (baw):  change norm on rayleigh part
c     00/10/27 (baw) created from stokespeel6
c     99/03/22 (baw) use rejection method only; clean up.
c     95/01/17 (mjw): add calls to ERRMSG,WRIMSG
c     
      
c     This subroutine
c     calculates scattering using analytic phase functions, either
c     Rayleigh scattering or Henyey-Greenstein.
C     The scattering matrix has four unique elements placed the following
c     way:
c     |  P1 P2 0   0  | 
c     |  P2 P1 0   0  |
c     |  0  0  P3 -P4 |
c     |  0  0  P4  P3 |
c     In this routine, I'm using the notation of White (1977).
c     note this has the sign reversed on P4, compared to 
c     the notation of Bohren&Huffman.
      
c     This matrix is rotated into and out of the scattering frame to
c     get stokes parameters in observer's frame
      
      use dust_mod
      implicit none
      
      character cmsgnm*70,cmsger*50
c     include 'stokes.txt'
      include 'random.txt'
 
      real*8 sip,sqp,sup,svp,cost,sint,pi,r2p,temp
     1     ,hit,htot
      
      real*8 peak,bmu,xran,ri1,ri3,cosi3,sini3,costp,sintp,phip,cosb2
     1     ,b,sinbt,sini2,bott,cosdph,cosi2,sin2i3,sin2i2,cos2i3,cos2i2
     1     ,sin2,cos2,sin2cos1,cos2sin1,p1,p2,p3,p4,a11,a12,a13,a21,a22
     1     ,a23,a24,a31,a32,a33,a34,a42,a43,a44,cosi1,sini1,sin2i1
     1     ,cos2i1,a,si,sq,su,sv,rprob
      
      integer idust,idust2
      
c     calculate peak here
c     modified Cornette & Shanks function---make sure peak is known
      if (g(idust2).gt.0.0) then
         temp=g2p1(idust2)-twog(idust2)
         peak=1.5*(onemg2(idust2))/(2.+g2(idust2))*2./temp/sqrt(temp)
      else
         temp=g2p1(idust2)+twog(idust2)
         peak=1.5*(onemg2(idust2))/(2.+g2(idust2))*2./temp/sqrt(temp)
      end if
c     H-G function
c     if (g(idust2).gt.0.0) then
c     temp=g2p1(idust2)-twog(idust2)
c     peak=onemg2(idust2)/temp/sqrt(temp)*2.
c     else
c     temp=g2p1(idust2)+twog(idust2)
c     peak=onemg2(idust2)/temp/sqrt(temp)*2.
c      end if
c     print*,'g,peak',g,peak
c     icount=0
      
      peak=peak*1.01

c     costp=cost
c     sintp=sint
c     phip=phi
      
c     5    continue
c     htot=htot+1.d0
c     xran=ran()         
c     bmu=1.d0-2.d0*xran
      cosb2=bmu**2
      b=cosb2-1.d0
      if(idust.eq.0) then
c     p1=1.d0+cosb2
c     p2=b*pl(idust2)
c     p3=(2.d0*bmu)
c     p4=0.d0 
c     baw, 2003/12/26
         p1=0.75*(1.d0+cosb2)
         p2=0.75*b*1.d0 
         p3=0.75*(2.d0*bmu)
         p4=0.d0 
         peak=1.501d0
      else 
         call dustmat(p1,p2,p3,p4,bmu,cosb2,pi,idust2)
      end if
      
      if(abs(bmu).gt.1.d0) then
         write(cmsger,'(a,f14.11)')'bmu eq.1',bmu
         call ERRMSG('WARNING','STOKES',cmsger)
         if(bmu.gt.1.d0) then
            bmu=1.d0
            cosb2=1.d0
            b=0.d0
         else
            bmu=-1.d0
            cosb2=1.d0
            b=0.d0
         end if
      end if
      
      sinbt=sqrt(1.d0-cosb2)
      
c     xran=ran()         
c     ri1=r2p*xran
      
c     **** ri1 gt pi ****
      
      if (ri1.gt.pi) then
         ri3=r2p-ri1
         cosi3=cos(ri3)
         sini3=sin(ri3)
         sin2i3=2.d0*sini3*cosi3
         cos2i3=2.d0*cosi3**2-1.d0
         a11=p1
         a12=p2*cos2i3
         a13=p2*sin2i3
         rprob=(a11*sip+a12*sqp+a13*sup)/sip
         if(rprob.gt.peak) then
            peak=rprob
            write(cmsgnm,'(f14.11,a)')peak,' = peak'
            call WRIMSG('STOKES',cmsgnm)
         end if
c     xran=ran()         
c     if(peak*xran.gt.rprob) go to 5
c     hit=hit+1
         a=rprob
         
         if(bmu.eq.1.d0) then
            write(cmsger,'(a,f14.11)')'bmu eq.1',bmu
            call ERRMSG('WARNING','STOKES',cmsger)
            go to 10
         else if(bmu.eq.-1.d0) then
            write(cmsger,'(a,f14.11)')'bmu eq.-1',bmu
            call ERRMSG('WARNING','STOKES',cmsger)
            sini3=1.d0
            sini2=1.d0
            cosi3=0.d0
            cosi2=0.d0
c     cost=-costp
c     sint=sintp
         else 
c     anything with p on the end is the primed angle, which is
c     the angle from whence it came. oops--except cosbwp.
            
c     cost=costp*bmu+sintp*sinbt*cosi3
c     write(6,*)'costp,bmu,sintp,sinbt,cosi3,cost'
c     write(6,*)costp,bmu,sintp,sinbt,cosi3,cost
c     everytime we do an arcsin or arccos we have to test to see
c     if roundoff error has caused cos or sin to be slightly
c     greater than 1.  then we have to set them to 1.
            if (abs(cost).lt.1.d0) then
c     sint=abs(sqrt(1.d0-cost**2))
               sini2=sini3*sintp/sint
c     if (sini2.lt.0) write(6,*)'sini2 lt 0',sini2
               bott=sint*sinbt
               cosi2=costp/bott-cost*bmu/bott
               if (abs(cosi2).gt.1.0001d0) then
                  write(cmsger,'(a,f14.11)')'cosi2 big',cosi2
                  call ERRMSG('WARNING','STOKES',cmsger)
               end if
            else
               sint=0.d0
               call ERRMSG('WARNING','STOKES','sint,sini2=0')
               sini2=0.0d0
               if (cost.ge.1.d0) cosi2=-1.d0
               if (cost.le.-1.d0) cosi2=1.d0
            end if
         end if
         
c     cosdph=-cosi2*cosi3+sini2*sini3*bmu
c     if (abs(cosdph).gt.1.d0) then
c     if(abs(cosdph).gt.1.0001d0) then
c     write(cmsger,'(a,f14.11)')'abs(cosdph).gt.1.001',cosdph
c     call ERRMSG('WARNING','STOKES',cmsger)
c     end if
c     if(cosdph.gt.1.d0) then
c     cosdph=1.0d0
c     else 
c     cosdph=-1.0d0
c     end if
c     end if
c     phi=phip+acos(cosdph)
c     if(phi.gt.r2p) phi=phi-r2p
c     if(phi.lt.0.d0) phi=phi+r2p
      
         sin2i2=2.d0*sini2*cosi2
         cos2i2=2.d0*cosi2**2-1.d0
         sin2=sin2i2*sin2i3
         cos2=cos2i2*cos2i3
         sin2cos1=sin2i2*cos2i3
         cos2sin1=cos2i2*sin2i3

         a21=p2*cos2i2
         a22=p1*cos2-p3*sin2
         a23=p1*cos2sin1+p3*sin2cos1
         a24=-p4*sin2i2
         a31=-p2*sin2i2
         a32=-p1*sin2cos1-p3*cos2sin1
         a33=-p1*sin2+p3*cos2
         a34=-p4*cos2i2
         a42=-p4*sin2i3
         a43=p4*cos2i3
         a44=p3
 
c     **** ri1 lt pi ****
      else   
 
         cosi1=cos(ri1)
         sini1=sin(ri1)
         sin2i1=2.d0*sini1*cosi1
         cos2i1=2.d0*cosi1**2-1.d0
         a11=p1
         a12=p2*cos2i1
         a13=-p2*sin2i1
         rprob=(a11*sip+a12*sqp+a13*sup)/sip
         if(rprob.gt.peak) then
            peak=rprob
            write(cmsgnm,'(f14.11,a)')peak,' = peak'
            call WRIMSG('STOKES',cmsgnm)
         end if
c     xran=ran()         
c     if(peak*xran.gt.rprob) go to 5
c     hit=hit+1
         a=rprob
         
         if(bmu.eq.1.d0) then
            write(cmsger,'(a,f14.11)')'bmu eq.1',bmu
            call ERRMSG('WARNING','STOKES',cmsger)
            go to 10
         else if(bmu.eq.-1.d0) then
            write(cmsger,'(a,f14.11)')'bmu eq.-1',bmu
            call ERRMSG('WARNING','STOKES',cmsger)
            sini1=1.d0
            sini2=1.d0
            cosi1=0.d0
            cosi2=0.d0
c     cost=-costp
c     sint=sintp
         else
         
c     anything with p on the end is the primed angle, which is
c     the angle from whence it came. oops--except cosbwp.
            
c     cost=costp*bmu+sintp*sinbt*cosi1
c     write(6,*)'costp,bmu,sintp,sinbt,cosi1,cost'
c     write(6,*)costp,bmu,sintp,sinbt,cosi1,cost
c     everytime we do an arcsin or arccos we have to test to see
c     if roundoff error has caused cos or sin to be slightly
c     greater than 1.  then we have to set them to 1.
            if (abs(cost).lt.1.d0) then
c     sint=abs(sqrt(1.d0-cost**2))
               sini2=sini1*sintp/sint
c     if (sini2.lt.0) write(6,*)'sini2 lt 0',sini2
               bott=sint*sinbt
               cosi2=costp/bott-cost*bmu/bott
               if (abs(cosi2).gt.1.0001d0) then
                  write(cmsger,'(a,f14.11)')'cosi2 big',cosi2
                  call ERRMSG('WARNING','STOKES',cmsger)
               end if
            else   
c     sint=0.d0
               call ERRMSG('WARNING','STOKES','sint,sini2=0')
               sini2=0.0d0
               if (cost.ge.1.d0) cosi2=-1.d0
               if (cost.le.-1.d0) cosi2=1.d0
            end if
         end if
         
c     cosdph=-cosi1*cosi2+sini1*sini2*bmu
c     if (abs(cosdph).gt.1.d0) then
c     if(abs(cosdph).gt.1.0001d0) then
c     write(cmsger,'(a,f14.11)')'abs(cosdph).gt.1.0001',cosdph
c     call ERRMSG('WARNING','STOKES',cmsger)
c     end if
c     if(cosdph.gt.1.d0) then
c     cosdph=1.0d0
c     else
c     cosdph=-1.0d0
c     end if
c     end if
         
c     phi=phip-acos(cosdph)
c     if(phi.gt.r2p) phi=phi-r2p
c     if(phi.lt.0.d0) phi=phi+r2p
         
         sin2i2=2.d0*sini2*cosi2
         cos2i2=2.d0*cosi2**2-1.d0
         sin2=sin2i2*sin2i1
         cos2=cos2i2*cos2i1
         sin2cos1=sin2i2*cos2i1
         cos2sin1=cos2i2*sin2i1
         a21=p2*cos2i2
         a22=p1*cos2-p3*sin2
         a23=-p1*cos2sin1-p3*sin2cos1
         a24=p4*sin2i2
         a31=p2*sin2i2
         a32=p1*sin2cos1+p3*cos2sin1
         a33=-p1*sin2+p3*cos2
         a34=-p4*cos2i2
         a42=p4*sin2i1
         a43=p4*cos2i1
         a44=p3
      end if
 
      si=(a11*sip+a12*sqp+a13*sup)
      sq=(a21*sip+a22*sqp+a23*sup+a24*svp)
      su=(a31*sip+a32*sqp+a33*sup+a34*svp)
      sv=(a42*sqp+a43*sup+a44*svp)
      
c     write(6,*)'si,sq,su',si,sq,su,nscat
      
      sip=si
      sqp=sq
      sup=su
      svp=sv
c     cosp=cos(phi)
c     sinp=sin(phi)
      
 10   continue
      return
      end
 
c     ***********************************************************************


