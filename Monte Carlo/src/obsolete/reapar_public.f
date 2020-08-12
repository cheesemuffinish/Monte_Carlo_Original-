      SUBROUTINE REAPAR(cflpar,iopar,atname,dustname)

C***********************************************************************
C Subroutine REAPAR handles the reading of input parameters from the
C ".par" file, as well as elementary processing with those input
C parameters to generate arrays.
c
C Original version created by P.J.Flatau, Colorado State Univ.
C Modified by B.T.Draine, Princeton Univ. Obs.
C ported from program DDSCAT by m. wolff
c
C History:
C 95/01/17 (mjw): stripped out of DDSCAT and ported to BLOB
C 99/12/02 (mjw): stripped out of BLOB and ported to MCTHERM
c                 parameter OCCULT defined in calling routine (SETUP)
c 00/03/19 (mjw): added some Vger stuff (Vger position in units of rmax)
c 00/06/30 (baw): redo some spot input
c 00/07/27 (mjw): remove extraneous commas that drove lf95 crazy
c                   (occurrences of "write(*,*),something")
c 03/08/13 (baw): get rid of parameters that confuse outside users
C
C***********************************************************************
C

      use tts_mod
      use grid_mod
      use opacin_mod
      use vger_mod
      use out_mod
      use filt_mod
      use messages
      use random
      implicit none

C     .. Scalar Arguments ..
      character cflpar*(*),atname*80,dustname(5)*80
      integer iopar,id


C     .. Local Scalars ..
      character cline*70,cmsgnm*300
      character rep*10
      character*3 cbub,chole,climb,cpeel,cpoly,cstream,czmin,clucy
     $     ,cwarn,coutillum
      character*6 cshape
     1   ,cwrite,cspot,cveeg,cnspot,cdiskacc,crminsub,calpha
     1   ,cplanckst,cdiskwarp
     
     
      character(len=50) :: par_dir
     
C     ..
C     .. Local Arrays ..
C     ..
C     .. External Subroutines ..
      external errmsg, wrimsg
C     ..
C     .. Intrinsic Functions ..
C     ..
C     .. Data statements ..
C     ..
c     .. Include statements ..

      include 'tabl.txt'
      include 'spot.txt'

      character npeel_c*3

C***********************************************************************

      open(unit=iopar,file=cflpar,status='old')

      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)

C***********************************************************************
C ... Preliminaries
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)

      READ(IOPAR,FMT=*,ERR=40)np
      WRITE(CMSGNM,FMT='(i11,a)')np
     &     ,' = NP = Number of photons from central star'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)npmin
      WRITE(CMSGNM,FMT='(i11,a)')npmin
     &     ,' = NPMIN =   Minimum number of photons for Lucy temperatur
     &e calculation'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) n_iter_min
      WRITE(CMSGNM,FMT='(i11,a)') n_iter_min
     &     ,' = NIMIN =   Minimum number of Lucy temperature iterations'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)cpeel
      if ( (cpeel.eq.'YES').OR.(cpeel.eq.'yes') ) then
         ipeel = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IPEEL=1, Using peeling-off procedure'
      else
         ipeel = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IPEEL=0, Using standard MC (no peeling-off)'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)
 
      READ(IOPAR,FMT=*,ERR=40)clucy
      if ( (clucy.eq.'YES').OR.(clucy.eq.'yes') ) then
         ilucy = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'ILUCY=1, Using Lucy T-correction'
      else
         ilucy = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'ILUCY=0, Using Bjorkman & Wood T-Correction'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

c      READ(IOPAR,FMT=*,ERR=40) accfrac
c      WRITE(CMSGNM,FMT='(1pe13.6,a)') accfrac
c     &     ,' = ACCFRAC = fraction of accretion luminosity'
c      CALL WRIMSG('REAPAR',CMSGNM)

c      READ(IOPAR,FMT=*,ERR=40) mdot
c      WRITE(CMSGNM,FMT='(1pe13.6,a)') mdot
c     &     ,' = MDOTD = Disk accretion rate'
c      CALL WRIMSG('REAPAR',CMSGNM)


c      READ(IOPAR,FMT=*,ERR=40)npout
c      WRITE(CMSGNM,FMT='(i11,a)')npout
c     &     ,' = NPOUT = Number of photons from outside star'
c      CALL WRIMSG('REAPAR',CMSGNM)
      npout=0

      READ(IOPAR,FMT=*,ERR=40)i1
      WRITE(CMSGNM,FMT='(i10,a)')i1
     &     ,' = I1 = random number seed'
      CALL WRIMSG('REAPAR',CMSGNM)
      i1=-(abs(i1))
      READ(IOPAR,FMT=*,ERR=40)iwrite
      WRITE(CMSGNM,FMT='(i9,a)')iwrite
     &     ,' = IWRITE = no. of photons between timing calls'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)cwarn
      if ( (cwarn.eq.'YES').OR.(cwarn.eq.'yes') ) then
         iwarn = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IWARN=1, Printing warning & error messages'
      else
         iwarn = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IWARN=0, Not printing warning/error messages'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)
      
      call set_warnings(iwarn)
      call set_errors(iwarn)
      call set_debug(iwarn)

c      READ(IOPAR,FMT=*,ERR=40)ctherm
c      if ( (ctherm.eq.'YES').OR.(ctherm.eq.'yes') ) then
c         itherm = 1
c         WRITE(CMSGNM,FMT='(a)')
c     &        ,'ITHERM=1, Calculating thermal emission + scattering'
c      else
c         itherm = 0
c         WRITE(CMSGNM,FMT='(a)')
c     &        ,'ITHERM=0, Calculating ONLY scattering'
c      end if
c      CALL WRIMSG('REAPAR',CMSGNM)



c      READ(IOPAR,FMT=*,ERR=40)cveeg
c      if ( (cveeg.eq.'YES').OR.(cveeg.eq.'yes') ) then
c         iveeg = 1
c         WRITE(CMSGNM,FMT='(a)')
c     &        'IVEEG=1, calculating radiation field at point inside'
c      else
c         iveeg = 0
c         WRITE(CMSGNM,FMT='(a)')
c     &        'IVEEG=0, not calculating rad. field at point inside'
c      end if
c      CALL WRIMSG('REAPAR',CMSGNM)
      iveeg=0

c      READ(IOPAR,FMT=*,ERR=40)cwrite
c      if ( (cwrite.eq.'YES').OR.(cwrite.eq.'yes') ) then
c         iwriteim = 1
c         WRITE(CMSGNM,FMT='(a)')
c     &        'IWRITEIM=1, writing out all images'
c      else
c         iwriteim = 0
c         WRITE(CMSGNM,FMT='(a)')
c     &        'IWRITEIM=0, only writing peeling-off images'
c      end if
c      CALL WRIMSG('REAPAR',CMSGNM)


  
      READ(IOPAR,FMT=*,ERR=40) par_dir
      print*,trim(par_dir)//' = parfile directory'
  
      filter_dir = par_dir

c      print*, dthresh,
c     $ ' = density threshold for choosing dustname(1) and (2)'

      READ(IOPAR,FMT=*,ERR=40)coutillum
      if ( (coutillum.eq.'YES').OR.(coutillum.eq.'yes') ) then
         iout = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IOUT=1, Outside illumination by ISRF'
      else
         iout = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IOUT=0, no external illumination '
      end if
      CALL WRIMSG('reapar',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)isrf_scl
      WRITE(CMSGNM,FMT='(1pe13.6,a)')isrf_scl
     &     ,' = ISRF_SCL = scale factor for ISRF'
      CALL WRIMSG('reapar',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)isrf_Av
      WRITE(CMSGNM,FMT='(1pe13.6,a)')isrf_Av
     &     ,' = ISRF_AV = extinction of ISRF'
      CALL WRIMSG('reapar',CMSGNM)

c******  PAH properties *******
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)

      READ(IOPAR,FMT=*,ERR=40)fpah
      WRITE(CMSGNM,FMT='(1pe13.6,a)')fpah
     &     ,' = FPAH = fraction of mass in PAHs'
      CALL WRIMSG('reapar',CMSGNM)
      
      READ(IOPAR,FMT=*,ERR=40) dustname(5)
      dustname(5) = trim(par_dir)//'/'//trim(dustname(5))

C***********************************************************************
C ... Central source properties
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      
      READ(IOPAR,FMT=*,ERR=40)climb
      if ( (climb.eq.'YES').OR.(climb.eq.'yes') ) then
        limb = 1
        WRITE(CMSGNM,FMT='(a)')
     &        'LIMB=1, Limb-darkening for central source'
      else
         limb = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'LIMB=0, NO limb-darkening for central source'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)cplanckst
      if ( (cplanckst.eq.'YES').OR.(cplanckst.eq.'yes') ) then
         iplanckst = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IPLANCKST=1, using planck function for stellar emission'
      else
         iplanckst = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IPLANCKST=0, using stellar atmosphere file for star'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)
       
      READ(IOPAR,FMT=*,ERR=40) atname
      atname = trim(par_dir)//'/'//trim(atname)
      print*,trim(atname)//' = atmosphere file'

      READ(IOPAR,FMT=*,ERR=40)rstar
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rstar
     &     ,' = RSTAR = Stellar radius in Solar radii'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)tstar
      WRITE(CMSGNM,FMT='(1pe13.6,a)')tstar
     &     ,' = TSTAR = Blackbody temperature of central star (K)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)massc
      WRITE(CMSGNM,FMT='(1pe13.6,a)')massc
     &     ,' = MASSC = Mass of central star (for TSC properties)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) cspot
      if ( (cspot.eq.'YES').OR.(cspot.eq.'yes') ) then
         ispot = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'ISPOT=1, emission from star+spot'
	   spotflag = 1   !require to run spotset first time
      else
         ispot = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'ISPOT=0, uniform emission from star'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

c      READ(IOPAR,FMT=*,ERR=40)tspot
c      WRITE(CMSGNM,FMT='(1pe13.6,a)')tspot
c     &     ,' = TSTAR = Blackbody temperature of spot (K)'
c      CALL WRIMSG('REAPAR',CMSGNM)
c     will calculate later

c      READ(IOPAR,FMT=*,ERR=40)thspot
c      WRITE(CMSGNM,FMT='(1pe13.6,a)')thspot
c     &     ,' = THSPOT = radius of spot in degrees'
c      CALL WRIMSG('REAPAR',CMSGNM)
c     replaced by fspot

c      READ(IOPAR,FMT=*,ERR=40)thspotin
c      WRITE(CMSGNM,FMT='(1pe13.6,a)')thspotin
c     &     ,' = THSPOTIN = inner radius of spot (for ring)'
c      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)spotlat
      WRITE(CMSGNM,FMT='(1pe13.6,a)')spotlat
     &     ,' = SPOTLAT = latitude of spot (degrees)'
      CALL WRIMSG('REAPAR',CMSGNM)

c      READ(IOPAR,FMT=*,ERR=40)spotlon
c      WRITE(CMSGNM,FMT='(1pe13.6,a)')spotlon
c     &     ,' = SPOTLON = longitude of spot (degrees)'
c      CALL WRIMSG('REAPAR',CMSGNM)
      spotlon=0.

      READ(IOPAR,FMT=*,ERR=40) cnspot
      if ( (cnspot.eq.'ONE').OR.(cnspot.eq.'one') ) then
         nspot = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'NSPOT=1, one spot'
      else
         nspot = 2
         WRITE(CMSGNM,FMT='(a)')
     &        'NSPOT=2, two spots'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

C***********************************************************************
C ... Disk properties
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      
      READ(IOPAR,FMT=*,ERR=40) dustname(1)
      dustname(1) = trim(par_dir)//'/'//trim(dustname(1))
      print*,trim(dustname(1))//' = dust file(1)'
      
      READ(IOPAR,FMT=*,ERR=40) dustname(2)
      dustname(2) = trim(par_dir)//'/'//trim(dustname(2))
      print*,trim(dustname(2))//' = dust file(2)'

      READ(IOPAR,FMT=*,ERR=40) dthresh
      
      READ(IOPAR,FMT=*,ERR=40)massd
      WRITE(CMSGNM,FMT='(1pe13.6,a)')massd
     &     ,' = MASSD = Disk mass in solar masses'
      CALL WRIMSG('REAPAR',CMSGNM)
      print*,'massd',massd
      READ(IOPAR,FMT=*,ERR=40)rmaxd
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmaxd
     &     ,' = RMAXD = Maximum disk radius'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) crminsub
      if ( (crminsub.eq.'YES').OR.(crminsub.eq.'yes') ) then
         irminsub = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IRMINSUB=1, scaling RMINE,RMIND to dust sub. radius'
      else
         irminsub = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IRMINSUB=0, scaling RMINE,RMIND to Rstar'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)rmind_in
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmind_in
     &     ,' = RMIND = Minimum disk radius'
      CALL WRIMSG('REAPAR',CMSGNM)
c      READ(IOPAR,FMT=*,ERR=40)rddust
c      WRITE(CMSGNM,FMT='(1pe13.6,a)')rddust
c     &     ,' = RDDUST = Minimum disk opacity radius in Rstar'
c      CALL WRIMSG('REAPAR',CMSGNM)
      rmind=rmind_in
      rddust=rmind_in

c     for grid version of code, leave this out
c      READ(IOPAR,FMT=*,ERR=40)rmind2
c      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmind2
c     &     ,' = RMIND2 = Minimum dust radius of full-scale density'
c      CALL WRIMSG('REAPAR',CMSGNM)
c      READ(IOPAR,FMT=*,ERR=40)fractd
c      WRITE(CMSGNM,FMT='(1pe13.6,a)')fractd
c     &     ,' = FRACTD = density scale factor between rmind and rmind2'
c      CALL WRIMSG('REAPAR',CMSGNM)
c     change this to RDDUST2 to have less confusion
      rmind2=rddust
      fractd=1.d0

      READ(IOPAR,FMT=*,ERR=40) czmin
      if ( (czmin.eq.'YES').OR.(czmin.eq.'yes') ) then
         izmin = 1
         WRITE(CMSGNM,FMT='(a)')
     &   'IZMIN=1, calculating disk scale height based on hseq at Rsub'
      else
         izmin = 0
         WRITE(CMSGNM,FMT='(a)')
     &    'IZMIN=0, using user input for ZMIN (scale height at Rstar)'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) (z1(id), id=1,ndg)
      WRITE(CMSGNM,FMT='(1pe13.6,a)') z1(1)
     &     ,' = ZMIN = Scale height of disk at Rstar'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40) (a(id), id=1,ndg)
      WRITE(CMSGNM,FMT='(1pe13.6,a)') a(1)
     &     ,' = A = Disk density exponent (~ r^(-a))'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40) (b(id), id=1,ndg)
      WRITE(CMSGNM,FMT='(1pe13.6,a)') b(1)
     &     ,' = B = Disk scale height exponent (~ r^(b))'
      CALL WRIMSG('REAPAR',CMSGNM)
      
c     I couldn't get it to read them in right, so set here
c      do id=2,ndg
c         z1(id)=z1(1)
c         b(id)=b(1)
c         a(id)=a(1)
c      end do
      
c      print*,'z1',z1
c      print*,'b',b
c      print*,'a',a
      
      READ(IOPAR,FMT=*,ERR=40) cdiskacc
      if ( (cdiskacc.eq.'YES').OR.(cdiskacc.eq.'yes') ) then
         idiskacc = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IDISKACC=1, Including disk accretion luminosity'
      else
         idiskacc = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IDISKACC=0, not including disk accretion luminosity'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) calpha
      if (((calpha.eq.'YES').OR.(calpha.eq.'yes'))
     #     .and.idiskacc.eq.1) then
         ialpha = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IALPHA=1, using alpha parameter to calculate disk Mdot'
      else
         ialpha = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IALPHA=0, inputting disk accretion rate'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      if (ialpha.eq.1) then

         READ(IOPAR,FMT=*,ERR=40) alphad
         WRITE(CMSGNM,FMT='(1pe13.6,a)') alphad
     &        ,' = ALPHAD = alpha disk parameter'
         CALL WRIMSG('REAPAR',CMSGNM)
         
      else

         READ(IOPAR,FMT=*,ERR=40) mdotdisk
         WRITE(CMSGNM,FMT='(1pe13.6,a)') mdotdisk
     &        ,' = MDOTDISK = disk accretion rate '
         CALL WRIMSG('REAPAR',CMSGNM)

      end if

      READ(IOPAR,FMT=*,ERR=40) rtrunc
      WRITE(CMSGNM,FMT='(1pe13.6,a)') rtrunc
     &     ,' = RTRUNC = magnetosphere co-rotation radius'
      CALL WRIMSG('REAPAR',CMSGNM)
      
      READ(IOPAR,FMT=*,ERR=40) fspot
      WRITE(CMSGNM,FMT='(1pe13.6,a)') fspot
     &     ,' = FSPOT = fractional area of spot'
      CALL WRIMSG('REAPAR',CMSGNM)
      if (fspot.eq.0.and.idiskacc.eq.1) then 
        print*,'Error: fspot=0 and CDISKACC=yes'
        print*,'0 area hotspot is a singularity'
        print*,'maybe you meant fspot=1? '
        print*,'stopping program'
        stop
      end if

      READ(IOPAR,FMT=*,ERR=40) cdiskwarp
      if ( (cdiskwarp.eq.'YES').OR.(cdiskwarp.eq.'yes') ) then
         idiskwarp = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IDISKWARP=1, Including disk warp'
      else
         idiskwarp = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IDISKWARP=0, not including disk warp'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) wsc
      WRITE(CMSGNM,FMT='(1pe13.6,a)') wsc
     &     ,' = WSC = scale for disk warp (1 = double scaleheight)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) wexp
      WRITE(CMSGNM,FMT='(1pe13.6,a)') wexp
     &     ,' = WEXP = exponent for disk warp (cos**wexp)'
      CALL WRIMSG('REAPAR',CMSGNM)

c     reset some variables if disk mass is 0
      if (massd.eq.0.d0) then
         izmin=1
         z1=1.d0
         idiskacc=0
         idiskwarp=0
         print*,'since disk mass is zero, setting:'
         print*,'CZMIN=YES, zscale=1.0, CDISKACC=NO '
      end if

C***********************************************************************
C ... Envelope properties
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      
      READ(IOPAR,FMT=*,ERR=40) dustname(3)
      dustname(3) = trim(par_dir)//'/'//trim(dustname(3))
      print*,trim(dustname(3))//' = dust file(3)'

      READ(IOPAR,FMT=*,ERR=40) dustname(4)
      dustname(4) = trim(par_dir)//'/'//trim(dustname(4))
      print*,trim(dustname(4)),' = dust file(4)'

      READ(IOPAR,FMT=*,ERR=40)rmax
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmax
     &     ,' = RMAX = Maximum envelope radius in AU'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rmine_in
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmine_in
     &     ,' = RMINE = Minimum envelope radius (rstar or rsub)'
      CALL WRIMSG('REAPAR',CMSGNM)
      rmine=rmine_in

c     for grid version of code, leave this out
      READ(IOPAR,FMT=*,ERR=40)rmine2
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmine2
     &     ,' = RMINE2 = Minimum disk radius of full-scale density'
      CALL WRIMSG('REAPAR',CMSGNM)
      
      READ(IOPAR,FMT=*,ERR=40)fracte
      WRITE(CMSGNM,FMT='(1pe13.6,a)')fracte
     &   ,' = FRACTE = density scale factor bet. rmine1 and rmine2'
      CALL WRIMSG('REAPAR',CMSGNM)
c      rmine2=rmine
c      fracte=1.d0

      READ(IOPAR,FMT=*,ERR=40)rate
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rate
     &     ,' = RATE = Mass infall rate for TSC env.'//
     &     '(M_O/yr)'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rc
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rc
     &     ,' = RC = centrifugal radius for TSC env.'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)chole
      READ(IOPAR,FMT=*,ERR=40)cshape
      READ(IOPAR,FMT=*,ERR=40)rchole      
      READ(IOPAR,FMT=*,ERR=40)thet1
c      READ(IOPAR,FMT=*,ERR=40)thet2
      READ(IOPAR,FMT=*,ERR=40)ex1
c      READ(IOPAR,FMT=*,ERR=40)ex2
      READ(IOPAR,FMT=*,ERR=40)z01
c      READ(IOPAR,FMT=*,ERR=40)z02
      READ(IOPAR,FMT=*,ERR=40)exf
      READ(IOPAR,FMT=*,ERR=40)rhoconst1
c      READ(IOPAR,FMT=*,ERR=40)rhoconst2
      READ(IOPAR,FMT=*,ERR=40)rhoamb

      thet2=thet1
      ex2=ex1
      z02=z01
c      if (thet1.gt.89.999) rhoconst1=rhoamb
      rhoconst2=rhoconst1
      print*,'rhoconst1,rhoconst2',rhoconst1,rhoconst2
      stop

      if ( (chole.eq.'YES').OR.(chole.eq.'yes') ) then
         ihole = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IHOLE=1, Carving cavity out'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(f6.1,a)')thet1
     &        ,' = THET1 = Opening angle of cavity wall'
c         CALL WRIMSG('REAPAR',CMSGNM)
c         WRITE(CMSGNM,FMT='(f6.1,a)')thet2
c     &        ,' = THET2 = Opening angle of outer cavity wall'
c         CALL WRIMSG('REAPAR',CMSGNM)

         WRITE(CMSGNM,FMT='(f7.2,a)')ex1
     &        ,' = EX1 = Cavity wall exponent'
         CALL WRIMSG('REAPAR',CMSGNM)

         thetmu0=thet2
         if (cshape.eq.'STREAM') then
            ipoly=0
            istream=1
            WRITE(CMSGNM,FMT='(a)')
     &           'cavit is streamline shape'
            WRITE(CMSGNM,FMT='(1pe11.4,a)')rchole
     &           ,' = RCHOLE = Streamline hole size in AU'
            CALL WRIMSG('REAPAR',CMSGNM)
            WRITE(CMSGNM,FMT='(f6.2,a)')thetmu0
     &           ,' = THETMU0 = Opening angle of streamline hole (deg)'
         else
            ipoly = 1
            istream=0
            WRITE(CMSGNM,FMT='(a)')
     &           'Polynomial-shaped cavity'
            WRITE(CMSGNM,FMT='(f7.2,a)')ex1
     &           ,' = EX1 =Cavity wall exponent'
            CALL WRIMSG('REAPAR',CMSGNM)
         end if

         WRITE(CMSGNM,FMT='(1pe11.4,a)')z01
     &        ,' = Z01 = Height of wall at w=0'
         CALL WRIMSG('REAPAR',CMSGNM)
c         WRITE(CMSGNM,FMT='(1pe11.4,a)')z02
c     &        ,' = Z02 = Height of outer wall at w=0'
c         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(f7.2,a)')exf
     &        ,' = EXF = exponent for cavity density power-law'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')rhoconst1
     &        ,' = RHOCONST1 = Coefficient for cavity density'
         CALL WRIMSG('REAPAR',CMSGNM)
c         WRITE(CMSGNM,FMT='(1pe11.4,a)')rhoconst2
c     &        ,' = RHOCONST2 = Coefficient for outer cavity density'
c         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')rhoamb
     &        ,' = RHOAMB = ambient density'
         CALL WRIMSG('REAPAR',CMSGNM)

c         if ( (cstream.eq.'YES').OR.(cstream.eq.'yes') ) then
c            istream = 1
c            WRITE(CMSGNM,FMT='(a)')
c     &           'ISTREAM=1, streamline hole'
c            CALL WRIMSG('REAPAR',CMSGNM)
c         
c         else
c            istream = 0
c            WRITE(CMSGNM,FMT='(a)')
c     &           'ISTREAM=0'
c         end if
      else
         ihole = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IHOLE=0, NO cavity'
         CALL WRIMSG('REAPAR',CMSGNM)
      end if

c      READ(IOPAR,FMT=*,ERR=40)cbub
c      READ(IOPAR,FMT=*,ERR=40)nbub
c      READ(IOPAR,FMT=*,ERR=40)zbub1
c      READ(IOPAR,FMT=*,ERR=40)zbub2
c      READ(IOPAR,FMT=*,ERR=40)buboa
c      if ( (cbub.eq.'YES').OR.(cbub.eq.'yes') ) then
c         ibub = 1
c         WRITE(CMSGNM,FMT='(a)')
c     &        'IBUB=1, including bubble'
c         CALL WRIMSG('REAPAR',CMSGNM)
c         WRITE(CMSGNM,FMT='(f7.2,a)')nbub
c     &        ,' = NBUB = Exponent describing bubble shape'
c         CALL WRIMSG('REAPAR',CMSGNM)
c         WRITE(CMSGNM,FMT='(1pe11.4,a)')zbub1
c     &        ,' = ZBUB1 = height of bubble inner wall'
c         CALL WRIMSG('REAPAR',CMSGNM)
c         WRITE(CMSGNM,FMT='(1pe11.4,a)')zbub2
c     &        ,' = ZBUB2 = height of bubble outer wall'
c         CALL WRIMSG('REAPAR',CMSGNM)
c         WRITE(CMSGNM,FMT='(1pe11.4,a)')buboa
c     &        ,' = BUBOA = Opening angle of hole in bubble'
c      else
c         ibub = 0
c         WRITE(CMSGNM,FMT='(a)')
c     &        'IBUB=0, NO bubble'
c      end if
c      CALL WRIMSG('REAPAR',CMSGNM)

      cbub='NO'
      nbub=4
      zbub1=250.
      zbub2=300.
      buboa=0.0

C***********************************************************************
C ... Vger stuff
c
c XVG0,YVG0,ZVG0 = location of VGER in units of RMAX
c XSRC0,YSRC0,ZSRC0 = location of outer illum. (RMAX)
c ISRC0,QSRC0,USRC0,VSRC0 = normalized stokes for incident outside illum.
c
c      READ(IOPAR,FMT=*,ERR=40)CLINE
c      CALL WRIMSG(' ',CLINE)
c      READ(IOPAR,FMT=*,ERR=40)xvg0,yvg0,zvg0
c      READ(IOPAR,FMT=*,ERR=40)xsrc0,ysrc0,zsrc0
c      READ(IOPAR,FMT=*,ERR=40)isrc0,qsrc0,usrc0,vsrc0
c      if (CVEEG.eq.'YES') then
c         write(CMSGNM,fmt='(3(1pe11.4,1x),a)') xvg0,yvg0,zvg0,
c     &        ' - VGER location in units of RMAX'
c         CALL WRIMSG('REAPAR',CMSGNM)
c         write(CMSGNM,fmt='(3(1pe11.4,1x),a)') xsrc0,ysrc0,zsrc0,
c     &        ' - src illum. in units of RMAX'
c         CALL WRIMSG('REAPAR',CMSGNM)
c         write(CMSGNM,fmt='(4(1pe11.4,1x),a)')isrc0,qsrc0,usrc0,
c     &        vsrc0,
c     &       ' - incident stokes'
c         CALL WRIMSG('REAPAR',CMSGNM)
c      end if

      xvg0=0.4
      yvg0=0.4
      zvg0=-0.2
      xsrc0=1.4
      ysrc0=0.0
      zsrc0=1.2
      isrc0=1.0
      qsrc0=0.0
      usrc0=0.0
      
      
C***********************************************************************
C ... output data
c

      CALL WRIMSG(' ',CLINE)
      READ(IOPAR,FMT=*,ERR=40)CLINE

      READ(IOPAR,FMT=*,ERR=40) rep
      output_imfilt = rep == 'YES' .or. rep == 'yes'
      if(output_imfilt) then
         WRITE(CMSGNM,FMT='(a)')
     & 'IMFILT=1, Outputing images convolved with filters'
      else
         WRITE(CMSGNM,FMT='(a)') 'IMFILT=0'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) rep
      output_imcube = rep == 'YES' .or. rep == 'yes'
      if(output_imcube) then
         WRITE(CMSGNM,FMT='(a)')
     & 'IMCUBE=1, Outputing multi-wavelength data cubes'
      else
         WRITE(CMSGNM,FMT='(a)') 'IMFILT=0'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) rep
      output_imsparse = rep == 'YES' .or. rep == 'yes'
      if(output_imsparse) then
         WRITE(CMSGNM,FMT='(a)')
     & 'IMSPARSE=1, Outputing multi-wavelength sparse matrice'
      else
         WRITE(CMSGNM,FMT='(a)') 'IMSPARSE=0'
      end if
      CALL WRIMSG('REAPAR',CMSGNM)
            
      READ(IOPAR,FMT=*,ERR=40) npeel
      WRITE(CMSGNM,FMT='(i11,a)') npeel
     &     ,' = NPEEL = Number of apertures for output'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) nap
      WRITE(CMSGNM,FMT='(i11,a)') nap
     &     ,' = NAP = Number of apertures for output'
      CALL WRIMSG('REAPAR',CMSGNM)
      
      allocate(aperture2(nap))

      READ(IOPAR,FMT=*,ERR=40)rmaxi
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmaxi
     &     ,' = RMAXI = Image half-size in AU'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)apmin
      WRITE(CMSGNM,FMT='((f8.1,1x),a)')apmin
     &     ,' = APMIN = radius of smallest aperture in AU'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)apmax
      WRITE(CMSGNM,FMT='((f8.1,1x),a)')apmax
     &     ,' = APMAX = radius of largest aperture in AU'
      CALL WRIMSG('REAPAR',CMSGNM)

      WRITE(npeel_c,'(I0)') npeel
      
      allocate(thete_arr(npeel),coste_arr(npeel),sinte_arr(npeel))

      READ(IOPAR,FMT=*,ERR=40) thete_arr
      WRITE(CMSGNM,FMT='('//trim(npeel_c)//'(f7.2,1x),a35)') thete_arr
     &,' = THETE = theta angle(s) (deg) of high S/N image(s)/SED(s)'
      CALL WRIMSG('REAPAR',CMSGNM)

      allocate(phie_arr(npeel),cospe_arr(npeel),sinpe_arr(npeel))
 
      READ(IOPAR,FMT=*,ERR=40) phie_arr
      WRITE(CMSGNM,FMT='('//trim(npeel_c)//'(f7.2,1x),a30)') phie_arr
     &,' = PHIE = phi angle(s) (deg) of high S/N image(s)/SED(s)'
      CALL WRIMSG('REAPAR',CMSGNM)

      phie=0.
      
      CLOSE (IOPAR)
      RETURN
   40 CONTINUE
      CALL ERRMSG('FATAL','REAPAR',' Error reading .par file')
      RETURN
      END




