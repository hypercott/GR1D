! -*-f90-*-
subroutine map_profile(lprofile_name)

  use GR1D_module
  implicit none

  character*(*) lprofile_name
  integer profile_zones
  
  real*8 buffer, dmass, dx
  integer i,ibuffer
  
  integer keytemp,keyerr,eosflag
  
  real*8, allocatable :: pradius(:), &
       pmass(:),prho(:),ptemp(:), &
       ppress(:),peps(:),pvel(:),&
       pye(:),pomega(:)

  real*8,allocatable,save :: pradius_new(:)
  
  real*8, parameter :: kboltz_cgs = 1.380662d-16
  real*8 :: rmax

  rmax = x1(n1) / length_gf
  
! read profile      
  open(666,file=trim(lprofile_name),status='unknown', & 
       form='formatted',action='read')
  read(666,*) profile_zones
  
  allocate(pradius(profile_zones),pmass(profile_zones))
  allocate(prho(profile_zones),ppress(profile_zones))
  allocate(ptemp(profile_zones))
  allocate(peps(profile_zones))
  allocate(pvel(profile_zones))
  allocate(pye(profile_zones))
  allocate(pomega(profile_zones))

! new short format:  
  if(profile_type.eq.1) then
     if(do_rotation) then
        do i=1,profile_zones
           read(666,*) ibuffer,pmass(i),pradius(i),&
                ptemp(i),prho(i),pvel(i),pye(i), &
                pomega(i)
        enddo
     else
        do i=1,profile_zones
           read(666,*) ibuffer,pmass(i),pradius(i),&
                ptemp(i),prho(i),pvel(i),pye(i), &
                buffer
        enddo
     endif
  else if(profile_type.eq.2) then
     ! In this case, we read in the pressure from the
     ! profile, then reset the temperature to have
     ! the same pressure as in the initial profile.
     write(6,*) "Profile Type 2!"
     write(6,*) "Resetting Temperature base on Profile Pressure"
     if(do_rotation) then
        do i=1,profile_zones
           read(666,*) ibuffer,pmass(i),pradius(i),&
                ptemp(i),prho(i),pvel(i),pye(i), &
                pomega(i),ppress(i)
        enddo
     else
        do i=1,profile_zones
           read(666,*) ibuffer,pmass(i),pradius(i),&
                ptemp(i),prho(i),pvel(i),pye(i), &
                buffer,ppress(i)
        enddo
     endif
     ! reset temperature based on pressure
     call map_profile_reset_temp(profile_zones,prho,ptemp,pye,ppress,&
          rmax,pradius)
  else
     write(6,*) "Profile type unknown!"
     stop
  endif

  ! go to c=G=Msun=1
  if (pmass(profile_zones).lt.1.0d10) then
  else
     pmass = pmass * mass_gf
  endif
  
  !is pradius the cell center or cell interface?  depends on input
  !models, GR1D assumes pradius is cell centers (along with the
  !density, etc.) unless you set WHW02profile to 1, then the code
  !below executes.  GR1D does not use the mass coordinate for creating
  !the initial profile

  if (WHW02) then
     allocate(pradius_new(profile_zones))
     pradius_new(1) = pradius(1)/2.0d0
     do i=2,profile_zones
        pradius_new(i) = (pradius(i)+pradius(i-1))/2.0d0
     enddo
     pradius(:) = pradius_new(:)
     !end of suggested code
  endif

  pradius = pradius * length_gf
  prho = prho * rho_gf
  pvel = pvel / clite
  pomega = pomega / time_gf
  
  if(eoskey.eq.3) then
     ptemp = ptemp * kboltz_cgs / mev_to_erg
  endif
  
  do i=1,n1
     call map_map(rho(i),x1(i),prho,pradius,profile_zones)
     call map_map(eps(i),x1(i),peps,pradius,profile_zones)
     call map_map(press(i),x1(i),ppress,pradius,profile_zones)
     call map_map(temp(i),x1(i),ptemp,pradius,profile_zones)
     call map_map(v1(i),x1(i),pvel,pradius,profile_zones)
     call map_map(ye(i),x1(i),pye,pradius,profile_zones)
     if (do_rotation) then
        call map_map(omega(i),x1(i),pomega,pradius,profile_zones)
     endif
  enddo

  if (do_rotation.and.set_omega) then
     !omega_c and omega_A already in code units
     omega(:) = omega_c/(1.0d0+(x1(:)/omega_A)**2)
  endif


  if(do_rotation) then
     write(*,*) "Have Rotation"
  ! set up vphi
     if(GR) then
        vphi(:) = omega(:)*x1(:)
     else
        vphi1(:) = omega(:)*x1(:)
     endif
  endif

  ! set up atmosphere
  call atmosphere_init
  
  do i=1,n1
         
! call eos
! first: reset eps
     keytemp = 1 ! coming in with temperature
     eosflag = 4 ! we want eps to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),eps(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        write(*,*) keyerr
        write(*,"(i3,i5,1P10E15.6)") eoskey,i,rho(i),temp(i),ye(i),eps(i)
        stop "problem in initial data: eos: eps"
     endif
     
     
! now get new pressures
     keytemp = 0 ! coming in with eps
     eosflag = 1 ! we want press to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),press(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        stop "problem in initial data: eos: press"
     endif

     if(eoskey.eq.3) then
! now get the entropy for analysis purposes
! only do this if we are using the hot nuclear eos!
        keytemp = 1 ! coming in with temperature
        eosflag = 8 ! we want entropy to be reset
        keyerr = 0
        call eos(i,rho(i),temp(i),ye(i),eps(i),entropy(i), &
             keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
        if(keyerr.ne.0) then
           stop "problem in initial data: eos: entropy"
        endif
     endif
! speed of sound
     keytemp = 1 ! coming in with temperature
     eosflag = 6 ! we want cs2 to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),cs2(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        stop "problem in initial data: eos: entropy"
     endif
  enddo

!  call output_single(press/press_gf,"presstest.xg")


end subroutine map_profile

! **************************************************************
subroutine map_linterp(x1,x2,y1,y2,x,y)

! perform linear interpolation      
  implicit none

  real*8 slope,x1,x2,y1,y2,x,y

  if (x2.lt.x1) then
     stop "Error in linterp!"
  endif

  slope = (y2 - y1) / (x2 - x1)

  y = slope*(x-x1) + y1
 
end subroutine  map_linterp

! ***************************************************************

subroutine map_find_index(zones,array,goal,upper_index,lower_index)
  
! bisection search
  implicit none
  
  integer zones,i
  real*8 array(*)
  real*8 goal
  integer middle_index,upper_index,lower_index

  lower_index = 1
  upper_index = zones
  
  do while ( (upper_index - lower_index) .gt. 1 )
     middle_index = (lower_index + upper_index) * 0.5
     if ( (goal .ge. array(lower_index)) &
          .and. (goal .le. array(middle_index)) ) then
        upper_index = middle_index
     else
        if ( (goal .ge. array(middle_index)) &
             .and. (goal .le. array(upper_index)) ) then
           lower_index = middle_index
        endif
     endif
  enddo
      
end subroutine map_find_index

! ******************************************************************

subroutine map_map(point_value,point_radius0,parray,pradius,zones)

  implicit none
  
  real*8 point_value, point_radius, point_radius0
  real*8 pradius(*), parray(*)
  integer zones
  integer upper_index, lower_index

  point_radius = abs(point_radius0)
  
  if (point_radius .ge. pradius(1) .and. & 
       point_radius .lt. pradius(zones) )  then
     
     call map_find_index(zones,pradius,point_radius, &
          upper_index,lower_index)
     
     call map_linterp( pradius(lower_index),pradius(upper_index), &
          parray(lower_index), parray(upper_index),  & 
          point_radius, point_value )

  else if (point_radius .lt. pradius(1)) then
     ! linear extrapolation
     call map_linterp(pradius(1),pradius(2), & 
          parray(1),parray(2),point_radius,point_value)

  else if (point_radius .gt. pradius(zones)) then
     ! linear extrapolation
     call map_linterp(pradius(zones-1),pradius(zones), & 
          parray(zones-1),parray(zones),point_radius,point_value)
  endif
  
  
end subroutine map_map

! ******************************************************************

subroutine map_limit(lprofile_name)
  use GR1D_module, only: grid_rmax,rho_cut
  implicit none
      
  character*(*) lprofile_name
  integer profile_zones
  
  real*8 buffer, dmass, dx
  integer i,ibuffer
  
  integer keytemp,keyerr,eosflag
  
  real*8, allocatable :: pradius(:), &
       pmass(:),prho(:),ptemp(:), &
       ppress(:),peps(:),pvel(:),&
       pye(:),pomega(:)
  
  
  real*8 :: kboltz_cgs = 1.380662d-16

! read profile      
  open(666,file=trim(lprofile_name),status='unknown', & 
       form='formatted',action='read')
  read(666,*) profile_zones
  
  allocate(pradius(profile_zones),pmass(profile_zones))
  allocate(prho(profile_zones),ppress(profile_zones))
  allocate(ptemp(profile_zones))
  allocate(peps(profile_zones))
  allocate(pvel(profile_zones))
  allocate(pye(profile_zones))
  allocate(pomega(profile_zones))
  
  do i=1,profile_zones
     read(666,*) ibuffer,pmass(i),pradius(i),&
          ptemp(i),prho(i),pvel(i),pye(i), &
          buffer
  enddo
  
  do i=1,profile_zones
     if (prho(i).gt.rho_cut) then
        grid_rmax = pradius(i)
     endif
  enddo
  
  write(*,*) "Maximum Radius (at rho=",rho_cut,") set to: ", grid_rmax
  
  close(666)
  
  deallocate(pradius,pmass)
  deallocate(prho,ppress)
  deallocate(ptemp)
  deallocate(peps)
  deallocate(pvel)
  deallocate(pye)
  deallocate(pomega)
  

end subroutine map_limit

! ******************************************************************

subroutine map_profile_isotopes

  use isomod
  use GR1D_module,only: n1,x1,ye,length_gf
  implicit none

  real*8, allocatable :: pradius(:)
  real*8, allocatable :: pcomp(:,:)

  integer :: profile_zones
  integer :: i,j
  real*8  :: normalize
  real*8  :: buffer

  write(6,*) "***** Setting up Isotope Advection"
  
  write(6,*) "Isotope file: ",trim(adjustl(isotope_profile))
  
  open(unit=666,file=trim(adjustl(isotope_profile)),form='formatted',&
       action='read')
  read(666,*) profile_zones, num_isotopes
  write(*,*) "We have ",profile_zones, "composition profile zones."
  write(*,*) "We have ",num_isotopes, "isotopes."

  allocate(pradius(profile_zones))
  allocate(pcomp(profile_zones,num_isotopes))
  
  ! allocate GR1D internal isotope datastructures
  call allocate_compositions(num_isotopes)

  !read in A's
  read(666,*) comp_details(1:num_isotopes,1)

  !read in Z's
  read(666,*) comp_details(1:num_isotopes,2)

  do i=1,profile_zones
     read(666,*) buffer,pradius(i),pcomp(i,1:num_isotopes)
  enddo
  close(666)

  !G = c = M_sun = 1
  pradius = pradius * length_gf

  ! map compositions onto grid
  do i=1,n1
     do j=1,num_isotopes
        call map_map(compositions(i,j),x1(i),pcomp(:,j),&
             pradius,profile_zones)
     enddo
  enddo
  
  ! normalize to one in each cell
  do i=1,n1
     normalize = sum(compositions(i,:))
     compositions(i,:) = compositions(i,:)/normalize
  enddo

  if(do_isotope_ye) then
     ! reset Ye based on isotope data
     ye(i) = sum(compositions(i,1:num_isotopes) * &
          (comp_details(1:num_isotopes,2) / &
           comp_details(1:num_isotopes,1) ))
     
  endif

  deallocate(pcomp)
  deallocate(pradius)
  write(6,*) "***** Done Setting up Isotope Advection"

end subroutine map_profile_isotopes

! ******************************************************************

subroutine map_profile_reset_temp(n,prho,ptemp,pye,ppress,&
     pradius,rmax)

  use GR1D_module,only: eoskey,eos_rf_prec,rho_gf,press_gf,&
       temp_mev_to_kelvin
  implicit none
  
  integer,intent(in)   :: n
  real*8,intent(in)    :: prho(*)
  real*8,intent(in)    :: pye(*)
  real*8,intent(in)    :: ppress(*)
  real*8,intent(in)    :: pradius(*)
  real*8,intent(in)    :: rmax
  real*8,intent(inout) :: ptemp(*)

  logical :: done
  integer :: i,eosflag,keyerr,keytemp
  integer :: count
  real*8  :: xrho,xtemp,xtemp0,xtemp1,xtemp2
  real*8  :: xpress,xpress0,xpress1,xpress2
  real*8  :: err,dummy,err1,err2
  real*8, parameter  :: prec = 1.0d-8
  integer, parameter :: countmax = 100
  real*8, parameter   :: rho_min = 1.5d3

  keytemp = 1 
  eosflag = 1 
  keyerr = 0

  do i=1,n
     if(pradius(i).ge.1.5d0*rmax) cycle
     if(prho(i).le.rho_min) cycle
     xrho = prho(i)*rho_gf
     xtemp0 = ptemp(i)/temp_mev_to_kelvin
     xpress0 = ppress(i)*press_gf

     call eos(i,xrho,xtemp0,pye(i),dummy,xpress, &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        stop "problem in initial data: map_profile_reset_temp"
     endif
     
     ! temperature always increases when P increases;
     ! let's hope that this here will bracket root
     err = xpress-xpress0
     if(err.lt.0) then
        xtemp1 = xtemp0
        xtemp2 = 2.0d0*xtemp0
        call eos(i,xrho,xtemp2,pye(i),dummy,xpress2, &
             keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
        if(keyerr.ne.0) then
           stop "problem in initial data: map_profile_reset_temp"
        endif

        if( (xpress-xpress0) * (xpress2-xpress0) .gt. 0.0d0) then
           write(6,*) "Branch 1"
           write(6,"(1P10E15.6)") xrho/rho_gf,xtemp0,xtemp2,pye(i),&
                (xpress2-xpress0)/xpress0, (xpress-xpress0)/xpress0
           stop "Not bracketing pressure root in map_profile!"
        endif

     else
        xtemp1 = 0.5d0*xtemp0
        xtemp2 = xtemp0
        call eos(i,xrho,xtemp1,pye(i),dummy,xpress1, &
             keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
        if(keyerr.ne.0) then
           stop "problem in initial data: map_profile_reset_temp"
        endif

        if( (xpress-xpress0) * (xpress1-xpress0) .gt. 0.0d0) then
           write(6,*) "Branch 2"
           write(6,"(1P10E15.6)") xrho/rho_gf,xtemp0,xtemp1,pye(i),&
                (xpress2-xpress0)/xpress0, (xpress-xpress0)/xpress0
           stop "Not bracketing pressure root in map_profile!"
        endif

     endif
     
     xtemp = 0.5d0 * (xtemp1 + xtemp2)

     count = 0
     done = .false.
     do while(.not. done .and. count.lt.countmax)
        count = count + 1
        
        call eos(i,xrho,xtemp,pye(i),dummy,xpress, &
             keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
        if(keyerr.ne.0) then
           stop "problem in initial data: map_profile_reset_temp"
        endif
        
        err = abs(xpress - xpress0)/xpress0
        if (err < prec) then
           done = .true.
        endif

        if(.not.done) then

           if( (xpress-xpress0) .gt. 0 ) then
              xtemp2 = xtemp
           else
              xtemp1 = xtemp
           endif

           xtemp = 0.5d0 * (xtemp1 + xtemp2)

        endif

     enddo

!     write(6,"(i6,1P10E15.6)") i,xtemp0,xtemp,err
     if (count.ge.countmax) then
        stop "bad -- no convergence in temperature reset"
     endif

!     call eos(i,xrho,xtemp,pye(i),dummy,xpress, &
!          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
!     write(6,"(i5,1P10E15.6)") i,xpress,xpress0,xtemp,xtemp0,pye(i)

     ptemp(i) = xtemp*temp_mev_to_kelvin

  enddo
!  stop "eh"


end subroutine map_profile_reset_temp
