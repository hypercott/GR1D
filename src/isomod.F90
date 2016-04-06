module isomod

  implicit none

  character(len=256) :: isotope_profile
  integer :: num_isotopes

  real*8,allocatable :: comp_details(:,:)  
  real*8,allocatable :: compositions(:,:)
  real*8,allocatable :: compositionsp(:,:)
  real*8,allocatable :: compositionsm(:,:)
  real*8,allocatable :: comp_q(:,:),comp_qm(:,:),comp_qp(:,:)
  real*8,allocatable :: comp_qold(:,:)
  real*8,allocatable :: comp_q_hat_old(:,:)
  real*8,allocatable :: comp_q_hat(:,:)
  real*8,allocatable :: comp_flux_diff(:,:)


contains
  subroutine allocate_compositions(num)

    use GR1D_module, only: n1
    implicit none
    integer :: num

    num_isotopes = num
    
    allocate(comp_details(num_isotopes,3))
    allocate(compositions(n1,num_isotopes))
    allocate(compositionsp(n1,num_isotopes))
    allocate(compositionsm(n1,num_isotopes))
    allocate(comp_q(n1,num_isotopes))
    allocate(comp_qold(n1,num_isotopes))
    allocate(comp_qp(n1,num_isotopes))
    allocate(comp_qm(n1,num_isotopes))
    allocate(comp_q_hat(n1,num_isotopes))
    allocate(comp_q_hat_old(n1,num_isotopes))
    allocate(comp_flux_diff(n1,num_isotopes))

    comp_details(:,:)   = 0.0d0
    compositions(:,:)   = 0.0d0
    compositionsm(:,:)  = 0.0d0
    compositionsp(:,:)  = 0.0d0
    comp_q(:,:)         = 0.0d0
    comp_qold(:,:)      = 0.0d0
    comp_qp(:,:)        = 0.0d0
    comp_qm(:,:)        = 0.0d0
    comp_q_hat(:,:)     = 0.0d0
    comp_q_hat_old(:,:) = 0.0d0
    comp_flux_diff(:,:) = 0.0d0

  end subroutine allocate_compositions

end module isomod
