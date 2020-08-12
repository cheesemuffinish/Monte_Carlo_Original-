subroutine initarr_it()

  use type_nested_image
  use tts_mod
  use grid_mod
  implicit none

  integer :: i,peelid
  real :: freq_min,freq_max

  print*,'in initarr_it'

  ! STANDARD

  ! SEDs - real*8, dimension(n_freq,n_mu,n_ph,n_apmax,n_output)
  si = 0.0d0
  sq = 0.0d0
  su = 0.0d0
  sv = 0.0d0

  ! SED uncertainties - real*8, dimension(n_freq,n_mu,n_ph,n_apmax,n_output)
  si2 = 0.0d0
  sq2 = 0.0d0
  su2 = 0.0d0
  sv2 = 0.0d0

  ! Number of photons in SEDs - real*8, dimension(n_freq,n_mu,n_ph,n_apmax,n_output)
  nums = 0.0d0
 
  ! aveinc = 0.0d0

  ! PEELOFF

  if(ipeel==1) then

    ! SEDs
    ! (n_freq,n_apmax,n_output)
  
   ti=0.0
   tq=0.0
   tu=0.0
   tv=0.0

    ! SED uncertainties
    ! (n_freq,n_apmax,n_output)
  
    ti2=0.0
    tq2=0.0
    tu2=0.0
    tv2=0.0
    
    numt=0.0

    ! Images - real*8, dimension(npeel,n_x,n_x,n_bands)

    image_b_i=0.0
    image_b_q=0.0
    image_b_u=0.0
    image_b_v=0.0
  
    ! Monochromatic:
  
    image_m_i_split=0.0
    image_m_i=0.0
    image_m_q=0.0
    image_m_u=0.0
    image_m_v=0.0

    ! Nested images
    
    do peelid=1,npeel
       freq_min = log10(real(numin)*2.9979e14/1.2398d0)
       freq_max = log10(real(numax)*2.9979e14/1.2398d0)
       call nested_image_setup(image_ms_i(peelid),100,1.e6,10.,nfreq,freq_min,freq_max,2)
       call nested_image_setup(image_ms_q(peelid),100,1.e6,10.,nfreq,freq_min,freq_max,2)
       call nested_image_setup(image_ms_u(peelid),100,1.e6,10.,nfreq,freq_min,freq_max,2)
       call nested_image_setup(image_ms_v(peelid),100,1.e6,10.,nfreq,freq_min,freq_max,2)
    end do

  end if

  ! PHYSICS

  tdust=tdust2
  nabs=0
  dtauabs=0.d0

end subroutine initarr_it
