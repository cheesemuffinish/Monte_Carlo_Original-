      subroutine atmosinit(atname)
      
      implicit none
      
      integer nHnu,nHnumax
      parameter(nHnumax=25000)
      real atmosnu,Hnuint
      dimension atmosnu(nHnumax),Hnuint(nHnumax)
      common /atmos_interp/atmosnu,Hnuint,nHnu


      CHARACTER atname*50

c-----------------------------------------------------------------------
c      write(*,*) 'Filename with the model atmosphere spectrum'
c      read(*,'(a50)') name
c      name='../parfiles/flx4000kg35.par'
      call partintgr(atname,atmosnu,Hnuint,nHnu)
      
      return
      end

c **********************************************************************
        SUBROUTINE partintgr(name,w,p,n)
c ----------------------------------------------------------------------
c   Reads in a two-column array containing wavelengths and fluxes
c   Calculates partial integrals at every wavelength point
c   Normalizes the integrals on the bolometric flux
c   Uses real function tr2 for integration
c-----------------------------------------------------------------------
c                                [Anatoly Miroshnichenko, Mar. 1998]
c ======================================================================       
        use dust_mod
        implicit none
      
        integer nHnu,nHnumax
        parameter(nHnumax=25000)

        CHARACTER name*50
        CHARACTER*80 header1,header2,header3
        INTEGER i,n
        REAL w(nHnumax),winv(nHnumax),f(nHnumax),chif(nHnumax),
     +       p(nHnumax),fbol,tr2
c-----------------------------------------------------------------------
        open(1,file=name,status='old')
        read(1,*) header1
        read(1,*) header2
        read(1,*) header3
        i=1
 1      read(1,*,end=2) w(i),winv(i),f(i)

	    w(i)=1.24*winv(i)
            call opacset(w(i))
            chif(i)=kappa*f(i)

            if(i.eq.1) then
                p(i)=0.
            else
	        p(i)=tr2(f,winv,i)
            end if

            i=i+1

        goto 1

 2      close(1)
        
	n=i-1
        
	fbol=p(n)
        do i=1,n
            p(i)=p(i)/fbol
        end do

	kappaf=tr2(chif,winv,i)/fbol

        write(*,*) 'kappaf=',kappaf        
	
	return
        end

c **********************************************************************
        REAL FUNCTION tr2(a,b,n)
c ----------------------------------------------------------------------
c   Numerical integration using the trapezium rule
c-----------------------------------------------------------------------
c                                [Anatoly Miroshnichenko, Feb. 1998]
c ======================================================================
        implicit none
        INTEGER i,n
        REAL a(n),b(n)
c-----------------------------------------------------------------------
        tr2=0.
        do 1 i=2,n
1       tr2=tr2+(a(i)+a(i-1))*(b(i-1)-b(i))/2.
        return
        end
c **********************************************************************
      
      
      real function atmosfreq()

      implicit none

      integer i0,i1,i
      real xi,f,P0,P1,P

      integer nHnu,nHnumax
      parameter(nHnumax=25000)
      real atmosnu,Hnuint
      dimension atmosnu(nHnumax),Hnuint(nHnumax)
      common /atmos_interp/atmosnu,Hnuint,nHnu

      real ran2
      integer iseed,idum
      common /ranparm/ iseed,idum

      
      xi=ran()

      i0=1
      i1=nHnu

      P0=  Hnuint(i0)
      if(xi.lt.P0) then
          atmosfreq=atmosnu(1)
	  return
      end if

      P1=  Hnuint(i1)
      if (P1 .gt. 1.d0) then
          write(*,*) 'atmosfreq: P1>1, P1 = ',P1
	  P1=1.
      end if
      if (xi.ge.P1) then
          atmosfreq=atmosnu(nHnu)
          return
      end if

      do while((i1-i0).gt.1)

          i=(i0+i1)/2
          P=Hnuint(i)

	  if (P .gt. 1.d0) then
              write(*,*)'atmosfreq: P>1; i,P = ',i,P
	      P=1.
	  end if

          if (P .lt. P0) then
	      write(*,*) 'atmosfreq: P<P0; P0,P,P1 = ',P0,P,P1
	      P0=P
	      i0=i
	  end if
          if (P .gt. P1) then
	      write(*,*) 'atmosfreq: P>P1; P0,P,P1 = ',P0,P,P1
	      P1=P
	      i1=1
	  end if
	  
          if (P.gt.xi) then
              i1=i
              P1=P
          else
              i0=i
              P0=P
          end if

      end do

      if (P1.le.P0) then
          write(*,*) 'atmosfreq: P1<=P0; P0, P1 = ',P0,P1
	  stop
      end if
      
c      f=(xi-P0)/(P1-P0)
      if (xi.lt.P0) then
          write(*,*) 'atmosfreq: xi<P0, xi,P0 = ',xi,P0
          f=0.
      else if (xi.gt.P1) then
          write(*,*) 'atmosfreq: xi>P1, xi,P1 = ',xi,P1
          f=1.
      else if (P0.le.0.d0) then
          write(*,*) 'atmosfreq: P0<=0, P0,i0,i1 = ',P0,i0,i1
	  f=0.5
      else if (P1.le.0.d0) then
          write(*,*) 'atmosfreq: P1<=0, P1 = ',P1
	  stop
      else if (xi.le.0.d0) then
          write(*,*) 'atmosfreq: xi<=0, xi = ',xi
	  stop
      else
          f=log(xi/P0)/log(P1/P0)
      end if
      if ((f.lt.0.d0).or.(f.ge.1.d0)) then
          write(*,*) 'atmosfreq: f = ',f,' out of bounds'
      end if

      atmosfreq= atmosnu(i0) * (atmosnu(i1)/atmosnu(i0))**f

      return
      end

