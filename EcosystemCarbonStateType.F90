module EcosystemCarbonStateType
  implicit none
  integer,parameter :: ifoliage  = 1
  integer,parameter :: iwoody    = 2
  integer,parameter :: imetlit   = 3
  integer,parameter :: istrlit   = 4
  integer,parameter :: imicsoil  = 5
  integer,parameter :: islowsoil = 6
  integer,parameter :: ipasssoil = 7
  type, public :: ecosystem_carbonstate_type
     real(kind=8),dimension(:),allocatable :: Cpool
     real(kind=8),dimension(:),allocatable :: Cpool_out

     contains
     procedure, public :: Init
  end type ecosystem_carbonstate_type

  contains
  subroutine Init(this,lu_init)
    use ParameterType, only : npools
    class(ecosystem_carbonstate_type) :: this
    integer,intent(in) :: lu_init
    integer i

    allocate(this%Cpool     (1:npools))
    allocate(this%Cpool_out (1:npools))
    read(lu_init,*)(this%Cpool(i),i=1,npools)

  end subroutine Init
end module EcosystemCarbonStateType

 
