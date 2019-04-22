module ParameterType
  
  implicit none
  integer,parameter :: npools=7
  type, public :: parameter_type
     real(kind=8) frac_foliage_to_metlit
     real(kind=8) frac_foliage_to_strlit
     real(kind=8) frac_woody_to_strlit
     real(kind=8) frac_metlit_to_micsoil
     real(kind=8) frac_strlit_to_micsoil
     real(kind=8) frac_slowsoil_to_micsoil
     real(kind=8) frac_passsoil_to_micsoil
     real(kind=8) frac_micsoil_to_slowsoil
     real(kind=8) frac_strlit_to_slowsoil
     real(kind=8) frac_micsoil_to_passsoil
     real(kind=8) frac_slowsoil_to_passsoil
     real(kind=8),dimension(1:npools) :: Calloc       
     real(kind=8),dimension(1:npools) :: k     
     contains
     procedure, public :: readparam
  end type parameter_type

  contains
  subroutine readparam(this,lu_param)
    class(parameter_type) :: this
    integer,intent(in)               :: lu_param
    integer i,j
    
    read(lu_param,*)(this%Calloc(i),i=1,npools)
    read(lu_param,*)(this%k(i),i=1,npools)
    read(lu_param,*)this%frac_foliage_to_metlit
    read(lu_param,*)this%frac_foliage_to_strlit
    read(lu_param,*)this%frac_woody_to_strlit
    read(lu_param,*)this%frac_metlit_to_micsoil
    read(lu_param,*)this%frac_strlit_to_micsoil
    read(lu_param,*)this%frac_slowsoil_to_micsoil
    read(lu_param,*)this%frac_passsoil_to_micsoil
    read(lu_param,*)this%frac_micsoil_to_slowsoil
    read(lu_param,*)this%frac_strlit_to_slowsoil
    read(lu_param,*)this%frac_micsoil_to_passsoil
    read(lu_param,*)this%frac_slowsoil_to_passsoil
    
   ! Unit from (day)-1 to (second)-1
    this%k(:) = this%k(:) / 3600 / 24 
  end subroutine readparam
end module ParameterType
