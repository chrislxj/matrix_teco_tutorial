module SimulationTimeType

  implicit none
  type, public :: simulation_time_type
     integer(kind=8),private      :: year
     integer(kind=8),private      :: doy
     integer(kind=8),private      :: startyr
     integer(kind=8),private      :: nyr_forc
     real(kind=8)   ,private      :: dt
     contains
     procedure, public :: init
     procedure, public :: tick_time
     procedure, public :: get_curr_time
     procedure, public :: is_endof_forcing
     procedure, public :: is_startof_forcing
     procedure, public :: is_endof_year
     procedure, public :: is_startof_year
     procedure, public :: get_dtime
     
  end type simulation_time_type

  contains
  subroutine init(this,year,startyr,nyr_forc,dt)
    class(simulation_time_type) :: this
    integer(kind=8),intent(in)          :: year
    integer(kind=8),intent(in)          :: startyr
    integer(kind=8),intent(in)          :: nyr_forc
    real(kind=8),intent(in)     :: dt
    
    this%year        = year
    this%doy         = 1
    this%startyr     = startyr
    this%nyr_forc    = nyr_forc
    this%dt          = dt
    
  end subroutine init

  subroutine get_curr_time(this, year, doy)
    class(simulation_time_type) :: this
    integer(kind=8),intent(out) :: year
    integer(kind=8),intent(out) :: doy
    year   = this%year
    doy    = this%doy
  end subroutine get_curr_time

  subroutine tick_time(this)
    class(simulation_time_type) :: this
    
    this%doy = this%doy + 1
    if(this%doy .gt. 365)then
       print*,'year',this%year
       this%year = this%year + 1
       this%doy = this%doy - 365
    end if
  end subroutine tick_time
 
  logical function is_endof_forcing(this)
    class(simulation_time_type) :: this
    is_endof_forcing = mod(this%year-this%startyr+1,this%nyr_forc) .eq. 0 .and. this%is_endof_year()
  end function is_endof_forcing
  
  logical function is_startof_forcing(this)
    class(simulation_time_type) :: this
    is_startof_forcing = mod(this%year-this%startyr,this%nyr_forc) .eq. 0 .and. this%is_startof_year()
  end function is_startof_forcing
    
  logical function is_endof_year(this)
    class(simulation_time_type) :: this
    is_endof_year = this%doy .eq. 365
  end function is_endof_year
    
  logical function is_startof_year(this)
    class(simulation_time_type) :: this
    is_startof_year = this%doy .eq. 1
  end function is_startof_year
    
  real(kind=8) function get_dtime(this)
    class(simulation_time_type) :: this
    get_dtime = this%dt
  end function get_dtime
  
end module SimulationTimeType
