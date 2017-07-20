!-*-f90-*-
!Initialize the M1 module
subroutine M1_init
  
  use GR1D_module, only : opacity_table,nulib_opacity_gf,nulib_emissivity_gf, &
       energy_gf, mev_to_erg,M1_maxradii,M1_imaxradii,length_gf,n1,x1,ghosts1,&
       ye,temp,rho_gf,rho_gf_inv,x1i,number_groups,rho,eas,ghosts1,q_M1,& 
       nulib_energy_gf, &
       temp_mev_to_kelvin,number_species,volume,pi,M1_moment_to_distro,clite, &
       hbarc_mevcm,M1_testcase_number,v_order,include_nes_kernels, &
       M1_moment_to_distro_inverse,nulib_kernel_gf,number_species_to_evolve, &
       include_epannihil_kernels,eos_rf_prec, sedonu, v1, vphi1
  use nulibtable

  implicit none
  include 'mpif.h'

  integer :: i, j, k, l
  real*8 :: tmp1, tmp2
  integer :: keytemp,keyerr
  real*8 :: eosdummy(17)
  real*8 :: xrho,xtemp,xye,xeta
  real*8 :: ensum
  integer myID, Nprocs, ierr
  
  call MPI_COMM_RANK (MPI_COMM_WORLD, myID, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, Nprocs, ierr)

  call nulibtable_reader(opacity_table,include_nes_kernels,&
       include_epannihil_kernels)
  
  !change units of emissivities, opacities, energies and inverse
  !energies to code units, then we only have to do it once these are
  !located in GR1D_module.F90

  !emissivities are in units of erg/cm^3/s/srad 
  !(I've multiplied by bin width (in MeV) in the table reader)
  !nulib_emissivity_gf = energy_gf/(length_gf**3*time_gf)

  !opacities are in units of cm^-1
  !nulib_opacity_gf = 1.0d0/length_gf

  !energies are in units of MeV
  !nulib_energy_gf = mev_to_erg*energy_gf

  !kernels are in units of cm^3 s^-1
  !nulib_kernel_gf = length_gf**3/time_gf

  nulibtable_emissivities = max(-200.0d0,& 
       log10(10.0d0**(nulibtable_emissivities)*nulib_emissivity_gf))
  nulibtable_absopacity = log10(10.0d0**nulibtable_absopacity*nulib_opacity_gf)
  nulibtable_scatopacity = log10(10.0d0**nulibtable_scatopacity* &
       nulib_opacity_gf)
  if (include_nes_kernels) then
     nulibtable_Itable_Phi0 = log10(10.0d0**nulibtable_Itable_Phi0* &
          nulib_kernel_gf)
     !Phi1 is stored as the ratio of Phi1/Phi0, so no log, no units
  endif
  if (include_epannihil_kernels) then
     nulibtable_epannihiltable_Phi0 = &
          log10(10.0d0**nulibtable_epannihiltable_Phi0*nulib_kernel_gf)
     !Phi1 is stored as the ratio of Phi1/Phi0, so no log, no units
  endif

  !convert energies to reduced units to save time, note nulib_ewidth
  !is NOT converted
  nulibtable_energies = nulibtable_energies*nulib_energy_gf
  nulibtable_ebottom = nulibtable_ebottom*nulib_energy_gf
  nulibtable_etop = nulibtable_etop*nulib_energy_gf
  nulibtable_inv_energies = nulibtable_inv_energies/nulib_energy_gf

  !set the log based version for time savings
  nulibtable_logenergies = log(nulibtable_energies)
  nulibtable_logetop = log(nulibtable_etop)

  !find zone that matches maximum radii
  do i=1,n1
     if (x1(i)/length_gf.lt.M1_maxradii) M1_imaxradii = i
  enddo

  ! check that density of outer M1 radius is above the minimum
  ! density of the NuLib table
  if (rho(M1_imaxradii)*rho_gf_inv.lt.&
       10.0d0**nulibtable_logrho(1) .and. myID==0) then
     write(6,*) "******************"
     write(6,*) "Density at maximum M1 radius is below NuLib table minimum!"
     write(6,"(A16,i5)") "M1_imaxradii = ",M1_imaxradii
     write(6,"(A24,1P1E15.6)") "radius(M1_imaxradii) = ", &
          x1(M1_imaxradii)/length_gf
     write(6,"(A21,1P1E15.6)") "rho(M1_imaxradii) = ",&
          rho(M1_imaxradii)*rho_gf_inv
     write(6,"(A17,1P1E15.6)") "min rho NuLib = ",&
          10.0d0**nulibtable_logrho(1)
     write(6,*) "Aborting!"
     write(6,*) "******************"
     stop
  endif

  if (v_order.eq.0.or.v_order.eq.-1) then
     if(myID==0) then
        write(*,*) "Velocity order is:",v_order," (-1 for all orders)"
     endif
  else
     stop "implement v_order"
  endif
  if(myID==0) then
     write(*,*) "M1_init: extract radii at", x1(M1_imaxradii)/length_gf, &
          "index:", M1_imaxradii, "of", n1
  endif

  !conversion from energy (momentum) density to angle integrated
  !distribution function (*\mu)
  M1_moment_to_distro(:) =  (2.0d0*pi*hbarc_mevcm)**3 / &
       (clite*nulib_emissivity_gf/nulib_opacity_gf * &
       nulibtable_ewidths(:)*mev_to_erg* &
       (nulibtable_energies(:)/nulib_energy_gf)**3)
  M1_moment_to_distro_inverse(:) = 1.0d0/M1_moment_to_distro(:)

  if (M1_imaxradii.gt.n1-ghosts1) then 
     if(myID==0) then
        write(6,*) "M1_init: Your extraction radii is too big"
        write(6,*) " "
     endif
     stop 
  endif

  if (number_species_to_evolve.eq.-1) then
     number_species_to_evolve = number_species
  endif

  call M1_updateeas

#ifdef HAVE_MC_CLOSURE
  call initialize_gr1d_sedonu(x1i/length_gf, n1,M1_imaxradii, ghosts1, sedonu)
#endif

  ! Note the following:
  ! In ZelmaniM1, we set radition initially to equilibrium.
  ! This makes sense at high density (i.e. in a proto-NS),
  ! where radiation and matter will actually be close to
  ! equilibrium. It doesn't make sense in a precollapse core.
  !
  ! By experiment, one finds that if radiation and matter are
  ! set to equilibrium initially, the nue energy density peaks
  ! at low densities and everything breaks in GR1D M1.
  !
  ! When turning on nua and nux later, it may make sense to
  ! set them to equilibrium at high density (but they'll go there
  ! quickly) and taper off to low density with a tanh taper
  ! function.

#if 0  
  ! The code below is by CDO and is debugging code
  ! that computes the equilibrium nue energy density and outputs
  ! it as a function of radius. Note that there seems to be something
  ! off with the units, I think.
  open(666,file="guender.dat")
  do i=ghosts1+1,M1_imaxradii
     xrho = rho(i) * rho_gf_inv
     xtemp = temp(i)
     xye = ye(i)
     keytemp = 1 !keep temperature                                                  
     keyerr = 0
     call nuc_eos_full(xrho,xtemp,xye,eosdummy(1),eosdummy(2),eosdummy(3), &
          eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),eosdummy(8), &
          eosdummy(9),eosdummy(10),eosdummy(11),eosdummy(12), &
          eosdummy(13),eosdummy(14),eosdummy(15),eosdummy(16),eosdummy(17), &
          keytemp,keyerr,eos_rf_prec)
     
     ! eosdummy(14) is mu_e
     ! eosdummy(17) is mu_hat = mu_n - mu_p
     ensum = 0.0d0
     xeta = (eosdummy(14)-eosdummy(17))/xtemp
     xeta = 0.0d0
     do j=1,number_groups
        tmp1 = nulibtable_energies(j)/(nulib_energy_gf*xtemp)
        tmp2 = clite*(nulibtable_energies(j)/nulib_energy_gf)**3 * &
             (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)/(1.0d0+exp(tmp1-xeta))
        ensum = ensum+tmp2*nulibtable_ewidths(j)
!        write(6,"(i5,1P10E15.6)") j,nulibtable_energies(j)/nulib_energy_gf,&
!             nulibtable_ewidths(j),ensum
     enddo
     if(myID==0) then
        write(6,"(i5,1P10E15.6)") i, x1(i)/length_gf, ensum*mev_to_erg
        write(666,"(i5,1P10E15.6)") i, x1(i)/length_gf, ensum*mev_to_erg
     endif

  enddo
  close(666)
#endif

end subroutine M1_init
