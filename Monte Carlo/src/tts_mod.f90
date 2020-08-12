module tts_mod

  use type_nested_image
  use grid_mod, only: ndg

  implicit none
  save

  character(len=10),parameter :: version = '20090821'

  ! STANDARD SEDS
  ! si/sq/su/sv     = used for fluxes
  ! si2/sq2/su2/sv2 = used for uncertainties
  ! nums            = number of photons used
  ! aveinc          = average inclination angle of the exiting photons in a theta bin--deleted

  integer :: nfreq ! number of frequencies for SEDs
  ! integer :: nfreqcube ! number of frequences for image cube
  ! 3-D outgoing angles
  ! integer,parameter :: nmu    = 20  ! number of mu viewing angles for SEDs
  ! integer,parameter :: nph    = 40  ! number of phi viewing angles for SEDs
  ! 2-D outgoing angles
  integer :: nmu ! number of mu viewing angles for SEDs
  integer :: nph ! number of phi viewing angles for SEDs
  ! 1=D outgoing angles
  ! integer,parameter :: nmu    = 2
  ! integer,parameter :: nph    = 1  ! number of phi viewing angles for SEDs
  !
  integer,parameter :: no     = 8   ! number of separate components
  ! (total/stellar/disk/envelope/outside/direct/scattered/thermal)
  integer,parameter :: no_init = 5 ! number of components of initial emitted luminosity sources
  ! (total/star not including hotspot/stellar hotspot/disk accretion/outside illumination)

  logical :: diffusion = .true.
  logical :: partial_peeloff = .false.

  ! (nfreq,nmu,nph,nap,no)
  real,allocatable,dimension(:,:,:,:,:) :: si,sq,su,sv,si2,sq2,su2,sv2,nums
  real,allocatable,dimension(:,:) :: si_init
  ! deleted average inclination array aveinc

  ! PEELOFF ANGLES

  integer :: npeel ! number of peel-off angles

  real*8,allocatable,dimension(:) :: thete_arr,phie_arr
  real*8,allocatable,dimension(:) :: coste_arr,sinte_arr
  real*8,allocatable,dimension(:) :: cospe_arr,sinpe_arr

  ! PEELOFF SEDS
  ! ti/tq/tu/tv     = used for fluxes
  ! ti2/tq2/tu2/tv2 = used for uncertainties
  ! numt            = number of photons used
  ! aveinc          = average inclination angle of the exiting photons in a theta bin--deleted
  ! nfreq,npeel,nap,no

  real,allocatable,dimension(:,:,:,:) :: ti,tq,tu,tv,ti2,tq2,tu2,tv2,numt

  ! PEELOFF IMAGES
  ! tihst/tqhst/tuhst/tvhst = used for fluxes

  ! nbnd needs to be equal to that in filt.txt
  integer,parameter :: nx    = 1
 ! integer,parameter :: nxhst = 149 ! image size
  integer :: nxhst ! image size
  integer,parameter :: nbnd  = 21  ! number of filters

  ! Broadband:
  ! (npeel,nxhst,nxhst,nbnd)

  real,allocatable,dimension(:,:,:,:) :: image_b_i,image_b_q,image_b_u,image_b_v

  ! Monochromatic:
  ! (npeel,nxhst,nxhst,nfreq,no)

  real,allocatable,dimension(:,:,:,:,:) :: image_m_i
  real,allocatable,dimension(:,:,:,:) :: image_m_q,image_m_u,image_m_v

  ! Sparse

  type(nested_image),allocatable,dimension(:) :: image_ms_i,image_ms_q,image_ms_u,image_ms_v

  real*8,allocatable :: tauenv(:,:,:)

  real*8,allocatable :: aperture2(:)
  real*8,allocatable :: u(:)
  real*8 :: b(ndg),a(ndg),rho0(ndg),fmass(ndg),taur(ndg),rhogap0(ndg)
  real*8 :: z1(ndg),agap(ndg),bgap(ndg),z1gap(ndg),rhoscale_gap(ndg)
  real*8 :: mdotcgs(ndg)
  real*8 :: nfrac(2)
  real*8 :: nscat,tot,zmax,rmax,lp,abflux
  real*8 :: rmin,rmind(ndg),rddust(ndg),rddustmin,flux,xmaxp
  real*8 :: sflux,aflux,cosb,sinb,rmaxd(ndg),rmaxdmax,rmaxi
  real*8 :: rpolbeam,dmu,dmuhalf,thetmu0,autors,rstar,wave
  real*8 :: tstar,tstarave,normstar,thete,phie,coste,sinte,cospe,sinpe
  real*8 :: sfactor,fractxh,numax,numin,lnurat,nunorm
  ! real*8 :: numaxcube,numincube,lnuratcube,nunormcube,nuratcube
  real*8 :: nurat,nub,accfrac,accfrac2,scount,dscount
  real*8 :: tshock,rtrunc,fspot,mdotdisk,apmax,apmin,ltot
  real*8 :: dthresh,l_isrf,isrf_Av,isrf_scl,tdiffavesave
  real*8 :: sigsave,isrf_frac,rmin_sg
  real*8 :: rhodens1,envexp
  real*8 :: rgapd1,rgapd2,rgape1,rgape2,fracte,rhogape
  real*8 :: pitch,sw,rspiral1,rspiral2 
  real*8 :: densratio,fractl,rimheight,gapheight,rimlength,gaplength,warpheight,warplength,radexp
  real*8 :: mdot,alphad,rimcurveheight,rimcurvelength,gapcurveheight,gapcurvelength,rimcurveexp
  real*8 :: gapcurveexp,radmisalign,incmisalign,lumout,sumall,sumsub

  integer :: np,npmin,n_iter_min,n_iter_max,npout,npbkg,iout
  integer :: npsav,limb,occult,nfinei,nfined
  integer :: nri,nzi,iwrite,isot,itherm,ipeel,ispot,i_orig,idiskacc
  integer :: ialpha,irminsub,izmin,iplanckst,nap,idiskwarp,ilucy
  integer :: ispiral,nspiral,ienvtype,ifractal,irimpuff,igappuff,iradexp
  integer :: igapd,igape,igapedens,SN,ihseq,iviscous,iter_hseq,nhseq1,nhseq2
  integer :: irimcurve,igapcurve,idiskcurve,imisalign,itemp,ihydro

  logical :: output_imfilt,output_imcube,output_imsparse
  
  integer :: sparse_dim, sparse_factor
  real :: sparse_rmin, sparse_rmax

  integer :: imave = 0
  integer :: imct = 0
  
  character(len=10) :: czmin

  character(len=50) :: par_dir

end module tts_mod

