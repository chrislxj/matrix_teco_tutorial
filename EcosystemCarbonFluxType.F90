module EcosystemCarbonFluxType
  implicit none
  type, public :: ecosystem_carbonflux_type
     real(kind=8) Cinput
     real(kind=8) Cinput_out
     real(kind=8) diag_I
     real(kind=8) diag_I_out
     real(kind=8) diag_ResidenceTime_out
     real(kind=8) diag_Capacity_out
     real(kind=8) diag_Potential_out
     real(kind=8),dimension(:,:),allocatable :: matrix_A
     real(kind=8),dimension(:,:),allocatable :: matrix_K
     real(kind=8),dimension(:)  ,allocatable :: matrix_B
     real(kind=8),dimension(:)  ,allocatable :: dC_dt
     
     real(kind=8),dimension(:)  ,allocatable :: diag_ResidenceTime
     real(kind=8),dimension(:)  ,allocatable :: diag_Capacity
     real(kind=8),dimension(:)  ,allocatable :: diag_Potential
     
     
     contains
     procedure, public :: Init
  end type ecosystem_carbonflux_type
  
  contains 
  subroutine Init(this)
     use ParameterType, only : npools
     class(ecosystem_carbonflux_type) :: this
     integer i
    
     allocate(this%matrix_A(1:npools,1:npools))
     allocate(this%matrix_K(1:npools,1:npools))
     allocate(this%matrix_B(1:npools))
     allocate(this%dC_dt   (1:npools))
     
     allocate(this%diag_ResidenceTime(1:npools))
     allocate(this%diag_Capacity     (1:npools))
     allocate(this%diag_Potential    (1:npools))

     this%matrix_A(:,:) = 0.
     this%matrix_K(:,:) = 0.
     this%matrix_B(:)   = 0.
     do i=1,npools
        this%matrix_A(i,i) = -1.
     end do

     this%Cinput_out = 0.

     this%diag_I_out             = -9999.
     this%diag_ResidenceTime_out = -9999.
     this%diag_Capacity_out      = -9999.
     this%diag_Potential_out     = -9999.
     
  end subroutine Init
     
end module EcosystemCarbonFluxType

 
