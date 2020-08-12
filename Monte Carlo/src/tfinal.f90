subroutine tfinal(nphot,iter,converge)

  use output_mod

  use tts_mod
  use grid_mod
  implicit none

  integer,intent(in) :: nphot
  integer,intent(in) :: iter
  character(len=3),intent(inout) :: converge

  character(len=4) :: suffix
  character(len=100) :: filename

  integer :: ir,it,ip,id,idust2,count,n
  real*8 :: dusttemp,tdiffave,told,tdiff2ave,tnew,tave,sig,rhosum,trhosum

  print*,'in tfinal'

  write(suffix,'("_",I3.3)') iter


  tdiffave=0.d0
  tdiff2ave=0.d0
  n=0

  do id=1,ndg
    
    if(.not.is_sg(id)) then
     idust2 = id

     do ir=1,nrg-1
        do it=1,ntg-1
           do ip=1,npg-1


              told=tdust2(ir,it,ip,id)
              if (.not.diffus(ir,it,ip)) then
                 !  idust2=dustarr(ir,it,ip)

                 tnew=dusttemp(ir,it,ip,dtauabs(ir,it,ip,id),tdust(ir,it,ip,id),ltot,nphot,idust2)
                 tdust2(ir,it,ip,id)=tnew
              else
                 tdust2(ir,it,ip,id)=tdust(ir,it,ip,id)
              end if

              !if (tdust2(ir,it,ip,id)*11605.d0lt.3.d0) tdust2(ir,it,ip,id)=3.d0/11605.d0

              !if (densarr(ir,it,ip,id).eq.0.d0) tdust2(ir,it,ip,id)=3.d0/11605.d0

              ! reset density to 0 if T>1600

   !      if (tdust2(ir,it,ip,id)*11605.d0.ge.1600.d0.and.densarr(ir,it,ip,id).gt.0.d0) then
   !   testing, 20111229, BAW
             if (tdust2(ir,it,ip,id)*11605.d0.ge.1600.d0) then
   !             print*,'dust sublimated',tdust2(ir,it,ip,id)*11605.
   !              tdust2(ir,it,ip,id) = 1600.d0 / 11605.d0
                  if(.not.diffus(ir,it,ip)) then
 !	                tdust2(ir,it,ip,id) = 0.5d0 / 11605.d0
 !                   densarr(ir,it,ip,id)=0.d0
 !                   massarr(ir,it,ip,id)=0.d0
   !  testing 20120128
 	                tdust2(ir,it,ip,id) = 1550.d0 / 11605.d0
                    densarr(ir,it,ip,id)=0.5*densarr(ir,it,ip,id)
                    massarr(ir,it,ip,id)=0.5*massarr(ir,it,ip,id)	               
                 end if
              end if

              if (densarr(ir,it,ip,id).eq.0.d0) then
	            tdust2(ir,it,ip,id)=0.1/11605.d0
	            massarr(ir,it,ip,id)=0.d0
	         !   diffus(ir,it,ip)=0   !sets entire diffusion region to 0!
	            nodust(ir,it,ip,id)=.true.
	          endif

              tnew=tdust2(ir,it,ip,id)

              if (tnew.ne.told.and.tnew>0.1/11605.d0.and.told>0.1/11605.d0) then      !not all cells have density in them
                 tdiffave=abs((tnew-told)/tnew)+tdiffave
                 tdiff2ave=((tnew-told)/tnew)**2.d0+tdiff2ave
                 n=n+1
              end if

           end do
        end do
     end do
     end if

  end do

  if (converge.eq.'no') then

     if (n.gt.0) then
     	tdiffave=tdiffave/dble(n)
     	tdiff2ave=tdiff2ave/dble(n)
     else
		tdiffave=0.d0
		tdiff2ave=0.d0
	endif

     ! compare tdiffave and tdiff2ave to saved values.  when these
     ! converge, converge='yes'

     print*,'tdiffave',tdiffave
     sig=sqrt(tdiff2ave-tdiffave**2.d0)
     print*,'old, new std dev',sigsave,sig

     print*,'fractional change in tdiffave ',(tdiffavesave-tdiffave)/tdiffave

     print*,'fractional change in std dev ',(sigsave-sig)/sig

     if (abs(sigsave-sig)/sig.lt.0.25d0) converge='yes'

     sigsave=sig
     tdiffavesave=tdiffave

     print*,'tfinal: converge = ',converge

  end if

  ! calculate dust destruction radius
  ! comment out when satisfied (yes I will regrid someday)

  it=(ntg+1)/2
  ir=4
  do id=1,ndg
    if(.not.is_sg(id)) then
      print*,'calculating dust destruction radius for id=',id
      do while (nodust(ir,it,1,id))
         ir=ir+1
    !     print*,tdust(ir,it,1,id)*11605.d0,tdust2(ir,it,1,id)*11605.d0,ir,it
      end do
      print*, 'dust destruction radius ',ir,rarr(ir-1)
    end if
  end do
  
  tdust = tdust   * 11605.d0
  tdust2 = tdust2 * 11605.d0

 ! call output_grid('tarr2'//suffix,tdust2)
    call output_grid('tarr' //suffix,tdust2)
  call output_grid('darr' //suffix,densarr)
  call output_grid('darrtot' //suffix,sum(densarr,dim=4))
  call output_grid('dtauabs'//suffix,dtauabs)
  ! call output_grid('tarr2',tdust2)
    call output_grid('tarr',tdust2)
  call output_grid('darr',densarr)
  call output_grid('darrtot',sum(densarr,dim=4))
  call output_grid('nabs',nabs)
  call output_grid('dtauabs',dtauabs)
!density-weight temperature
  do ir=1,nrg-1
  do it=1,ntg-1
  do ip=1,npg-1
	rhosum=0.d0
	trhosum=0.d0
    do id=1,ndg
	  if (tdust2(ir,it,ip,id).gt.0.1/11605.d0.and.densarr(ir,it,ip,id).gt.0.d0) then
		rhosum=rhosum+densarr(ir,it,ip,id)
		trhosum=trhosum+densarr(ir,it,ip,id)*tdust2(ir,it,ip,id)
	  endif
    enddo
    tdust2dw(ir,it,ip)=trhosum/rhosum
  enddo
  enddo
  enddo  
  call output_grid('tarrdw' //suffix,tdust2dw)
  call output_grid('tarrdw',tdust2dw)

  tdust = tdust   / 11605.d0
  tdust2 = tdust2 / 11605.d0

  ! see how things are looking
  ! see how things are looking
  do id=1,ndg
     write(filename,'("tave1_",I0,".dat")') id
     open(unit=15,file=filename,status='unknown')
     write(15,*) 'index,r(rstar),r(au),Tave(r)'
     do ir=1,nrg-1
        tave=0.d0
        count=0
        do it=1,ntg-1
           do ip=1,npg-1
              tave=tdust2(ir,it,ip,id)+tave
              count=count+1
           end do
        end do
        tave=tave/dble(count)
        write(15,'(i10,3(2x,e13.6),3(f10.4))') ir,ravearr(ir),ravearr(ir)/autors,tave*11605.d0
     end do
     close(15)
  end do

  print*,'done with tfinal'

  return
end subroutine tfinal
