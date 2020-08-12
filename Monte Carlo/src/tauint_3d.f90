module tauint_mod
  implicit none
contains

  subroutine tauint(iphot,tsum,ir,it,ip)

    use tts_mod
    use grid_mod
    use taunum_mod
    use opacin_mod
    use messages
    use dust_mod
    use closest_wall, only : reset_t,find_next_wall
    use constants
    implicit none

    character(len=20),parameter :: routine = 'tauint_3d'

    real*8 :: tsum,dt,opac,dtau,ds
    integer :: iface,iphot,it,ir,ip,id

    ! set optical depth sum to zero
    tsum=0.d0

    do

       call testgrid(xp,yp,zp,rsq,rtot,ir,it,ip,on_wall,iphot,routine)

       call reset_t()

       ! find shortest non-negative distance to two constant radius surfaces
       call radface(ux,uy,uz,xp,yp,zp,r2arr(ir),r2arr(ir+1))

       ! find shortest positive distance to two constant theta surfaces
       if (ntg.gt.2) call thetaface(ux,uy,uz,xp,yp,zp,tan2arr(it),tan2arr(it+1))

       ! find shortest non-negative distance to two constant phi surfaces
       if (npg.gt.2.and.ntg.gt.2) call phiface(aarr(ip),barr(ip),0.d0,0.d0,aarr(ip+1),barr(ip+1),0.d0,0.d0,ux,uy,uz,xp,yp,zp)

       call find_next_wall(on_wall,dt,iface)

       if (dt.lt.0.d0) then
          print*,'problem with find_wall! did not find surface! '
          print*,'in tauint,ir,rtot2,r2arr(ir)',ir,rsq,r2arr(ir)
          print*,'ux,uy,uz,iphot',ux,uy,uz,iphot
          print*,'xtot,ytot,ztot',xp,yp,zp
          stop
       end if

       opac=0.d0

       do id=1,ndg
          ! in lucy method, densarr gets set to zero if T>Tsub so don't have to check this
          opac=opac+kapd(id)*densarr(ir,it,ip,id)
       end do

       dtau=opac*dt

       if ((tsum+dtau).gt.tau) then
          ! sum up path for mean intensity
          ds=(tau-tsum)/opac
          ! print *,ds,rsol,rstar
          ! print *,ds*rsol*rstar
 !         if (nub.gt.1.2398d0/0.55) then !I think this is lambda<0.55 microns
          ! print *,jmean
!          jmean(ir,it,ip)=jmean(ir,it,ip)+ds*rsol*rstar*kappa(5)*kappav(5)*(1.d0-albedo(5))
          jmean(ir,it,ip)=jmean(ir,it,ip)+ds*kapd(5)*(1.d0-albedo(5))
        !  jmean(ir,it,ip)=jmean(ir,it,ip)+ds*rsol*rstar
          ! dtauabs(ir,it,ip)=dtauabs(ir,it,ip)+ds*opac*(1.d0-albedo(idust2))
!        end if
          do id=1,ndg
             dtauabs(ir,it,ip,id)=dtauabs(ir,it,ip,id)+ds*kapd(id)*densarr(ir,it,ip,id)*(1.d0-albedo(id))
          end do
          xp=xp+ds*ux
          yp=yp+ds*uy
          zp=zp+ds*uz
          rsq=xp**2+yp**2+zp**2
          rtot=sqrt(rsq)
          on_wall = .false.
          ! we are done with optical depth integration
          return
       end if

       ! sum up path for mean intensity
!        if (nub.gt.1.2398d0/0.55) then !I think this is lambda<0.55 microns
!   jmean(ir,it,ip)=jmean(ir,it,ip)+dt*rsol*rstar*kappa(5)*kappav(5)*(1.d0-albedo(5))
      jmean(ir,it,ip)=jmean(ir,it,ip)+dt*kapd(5)*(1.d0-albedo(5))
     !  jmean(ir,it,ip)=jmean(ir,it,ip)+dt*rsol*rstar
       ! print*,'wave,nub',wave,nub,1.2398d0/nub
!        end if
       do id=1,ndg
          dtauabs(ir,it,ip,id)=dtauabs(ir,it,ip,id)+dt*kapd(id)*densarr(ir,it,ip,id)*(1.d0-albedo(id))
          ! if (id.eq.2) print*,'dtauabs ',dtauabs(ir,it,ip,id)
       end do
       ! dtauabs(ir,it,ip)=dtauabs(ir,it,ip)+dt*opac*(1.d0-albedo(idust2))

       select case(iface)
       case(1)
          ir=ir-1
       case(2)
          ir=ir+1
       case(3)
          it=it-1
       case(4)
          it=it+1
       case(5)
          ip=ip-1
       case(6)
          ip=ip+1
       case default
          stop "ERROR: invalid iface"
       end select

       ! account for crossing the phi=0 plane:
       if (ip == 0)   ip = npg - 1
       if (ip == npg) ip = 1

       xp=xp+dt*ux
       yp=yp+dt*uy
       zp=zp+dt*uz
       on_wall = .true.
       rsq=xp*xp+yp*yp+zp*zp
       rtot=sqrt(rsq)
       tsum=tsum+dtau

  !     print*, 'tsum,dtau,opac,dens',tsum,dtau,opac,densarr(ir,it,ip,3)
       
       ! if photon exists grid, exit integration
       if (ir == nrg) then
          exitflag=1
          return
       end if

       ! if photon intersects star, re-emit
       if (ir == 0) then
          aflux=aflux+1.d0
          rtot=sqrt(rsq)
          ! reemit thermal photon
          call initpabs
          call opacset(nub)
          wave=1.2398d0/nub
          i_orig=1
          iflag=0
          ! idust2=dustarr(1,1,1)
          do id=1,ndg
             kapd(id)=kappa(id)*rsol*rstar*kappav(id)
          end do
          ! rlam(idust2)=albedo(idust2)
          ir=1
          tsum=0.d0
          !set new random optical depth! 20090720,  JB and BAW
       end if

    end do

  end subroutine tauint

end module tauint_mod
