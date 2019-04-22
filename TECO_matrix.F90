Program CLM_TECOMATRIX

  use EcosystemCarbonStateType
  use EcosystemCarbonFluxType
  use ParameterType
  use SimulationTimeType
  implicit none

  integer,parameter :: lu_Cinput=100,lu_param=110,lu_init=120,lu_out=130,lnamelist=10
  integer(kind=8) iyr,startyr,year,doy,nyr_sim,nyr_forc
  real(kind=8) dt
  character(len=150) :: file_Cinput,file_param,file_init,file_out
  logical use_matrix,is_SASU_spinup,is_diagnostic
  logical is_write_yearly

  type(parameter_type)      :: param
  type(ecosystem_carbonstate_type) :: cs_ecosys
  type(ecosystem_carbonflux_type)  :: cf_ecosys
  type(simulation_time_type)       :: SimTime

  namelist /control_pm/ file_Cinput,file_param,file_init,file_out,nyr_sim,nyr_forc,&
                        startyr,dt,use_matrix,is_SASU_spinup,is_diagnostic,is_write_yearly
                        
  open(lnamelist,file='teco.nml')
  
  call system('cat teco.nml')
  
  read(lnamelist,nml=control_pm)

  call SimTime%init(startyr,startyr,nyr_forc,dt)

  call open_io_file(lu_Cinput,lu_param,lu_init,lu_out,file_Cinput,file_param,file_init,file_out)

  call param%readparam(lu_param)

  call cs_ecosys%init(lu_init)

  call cf_ecosys%init()

  do iyr = 1, nyr_sim*365
      
     call SimTime%get_curr_time(year,doy)
     
     if(SimTime%is_startof_forcing())then
        rewind(lu_Cinput)
     end if

     call C_input(lu_Cinput,cf_ecosys)

     call update_ecosysC(param,cf_ecosys,cs_ecosys,use_matrix,SimTime)
     
     if(is_diagnostic)then
         call matrix_diagnostic(cf_ecosys,cs_ecosys,SimTime)
     end if
     
     if(is_SASU_spinup)then
        call matrix_spinup(cf_ecosys,cs_ecosys,SimTime)
     end if
     
     call write_output(lu_out,cf_ecosys,cs_ecosys,SimTime,is_write_yearly,is_diagnostic)
     
     call SimTime%tick_time()
  end do

end Program CLM_TECOMATRIX

subroutine update_ecosysC(param,cf_ecosys,cs_ecosys,use_matrix,SimTime)
  use EcosystemCarbonStateType, only : ecosystem_carbonstate_type,&
          ifoliage,iwoody,imetlit,istrlit,imicsoil,islowsoil,ipasssoil
  use EcosystemCarbonFluxType,  only : ecosystem_carbonflux_type
  use ParameterType, only : parameter_type,npools
  use SimulationTimeType, only : simulation_time_type

  type(parameter_type)            ,intent(in)    :: param
  type(ecosystem_carbonflux_type) ,intent(inout) :: cf_ecosys
  type(ecosystem_carbonstate_type),intent(inout) :: cs_ecosys
  logical                         ,intent(in)    :: use_matrix
  type(simulation_time_type)      ,intent(in)    :: SimTime
  

!!local variables

  real(kind=8) dt

  dt = SimTime%get_dtime()
  
! plant
  if(.not. use_matrix)then
     cf_ecosys%dC_dt(ifoliage)  = cf_ecosys%Cinput * param%Calloc(ifoliage) &
                                - param%k(ifoliage) * cs_ecosys%Cpool(ifoliage)
     cf_ecosys%dC_dt(iwoody)    = cf_ecosys%Cinput * param%Calloc(iwoody)   &
                                - param%k(iwoody)   * cs_ecosys%Cpool(iwoody)

! litter
     cf_ecosys%dC_dt(imetlit)   = cs_ecosys%Cpool(ifoliage)  * param%k(ifoliage) * param%frac_foliage_to_metlit &
                                - cs_ecosys%Cpool(imetlit)   * param%k(imetlit)  
     cf_ecosys%dC_dt(istrlit)   = cs_ecosys%Cpool(ifoliage)  * param%k(ifoliage) * param%frac_foliage_to_strlit &
                                + cs_ecosys%Cpool(iwoody)    * param%k(iwoody)   * param%frac_woody_to_strlit   &
                                - cs_ecosys%Cpool(istrlit)   * param%k(istrlit)

! soil
     cf_ecosys%dC_dt(imicsoil)  = cs_ecosys%Cpool(imetlit)   * param%k(imetlit)   * param%frac_metlit_to_micsoil    &
                                + cs_ecosys%Cpool(istrlit)   * param%k(istrlit)   * param%frac_strlit_to_micsoil    &
                                + cs_ecosys%Cpool(islowsoil) * param%k(islowsoil) * param%frac_slowsoil_to_micsoil  &
                                + cs_ecosys%Cpool(ipasssoil) * param%k(ipasssoil) * param%frac_passsoil_to_micsoil  &
                                - cs_ecosys%Cpool(imicsoil)  * param%k(imicsoil)
     cf_ecosys%dC_dt(islowsoil) = cs_ecosys%Cpool(imicsoil)  * param%k(imicsoil)  * param%frac_micsoil_to_slowsoil  &
                                + cs_ecosys%Cpool(istrlit)   * param%k(istrlit)   * param%frac_strlit_to_slowsoil   &
                                - cs_ecosys%Cpool(islowsoil) * param%k(islowsoil) 
     cf_ecosys%dC_dt(ipasssoil) = cs_ecosys%Cpool(imicsoil)  * param%k(imicsoil)  * param%frac_micsoil_to_passsoil  &
                                + cs_ecosys%Cpool(islowsoil) * param%k(islowsoil) * param%frac_slowsoil_to_passsoil &
                                - cs_ecosys%Cpool(ipasssoil) * param%k(ipasssoil)  
     
     cs_ecosys%Cpool(:)  = cs_ecosys%Cpool(:)  + cf_ecosys%dC_dt(:) * dt

  else
     cf_ecosys%matrix_B(ifoliage) = param%Calloc(ifoliage)
     cf_ecosys%matrix_B(iwoody)   = param%Calloc(iwoody)

     cf_ecosys%matrix_K(ifoliage ,ifoliage)  = param%k(ifoliage )
     cf_ecosys%matrix_K(iwoody   ,iwoody  )  = param%k(iwoody   )
     cf_ecosys%matrix_K(imetlit  ,imetlit )  = param%k(imetlit  )
     cf_ecosys%matrix_K(istrlit  ,istrlit )  = param%k(istrlit  )
     cf_ecosys%matrix_K(imicsoil ,imicsoil)  = param%k(imicsoil )
     cf_ecosys%matrix_K(islowsoil,islowsoil) = param%k(islowsoil)
     cf_ecosys%matrix_K(ipasssoil,ipasssoil) = param%k(ipasssoil)

     cf_ecosys%matrix_A(imetlit  ,ifoliage)  = param%frac_foliage_to_metlit
     cf_ecosys%matrix_A(istrlit  ,ifoliage)  = param%frac_foliage_to_strlit
     cf_ecosys%matrix_A(istrlit  ,iwoody)    = param%frac_woody_to_strlit
     cf_ecosys%matrix_A(imicsoil ,imetlit)   = param%frac_metlit_to_micsoil
     cf_ecosys%matrix_A(imicsoil ,istrlit)   = param%frac_strlit_to_micsoil
     cf_ecosys%matrix_A(imicsoil ,islowsoil) = param%frac_slowsoil_to_micsoil
     cf_ecosys%matrix_A(imicsoil ,ipasssoil) = param%frac_passsoil_to_micsoil
     cf_ecosys%matrix_A(islowsoil,imicsoil)  = param%frac_micsoil_to_slowsoil
     cf_ecosys%matrix_A(islowsoil,istrlit)   = param%frac_strlit_to_slowsoil
     cf_ecosys%matrix_A(ipasssoil,imicsoil)  = param%frac_micsoil_to_passsoil
     cf_ecosys%matrix_A(ipasssoil,islowsoil) = param%frac_slowsoil_to_passsoil

     cf_ecosys%dC_dt(:) = cf_ecosys%Cinput * cf_ecosys%matrix_B(:) &
                        + matmul(matmul(cf_ecosys%matrix_A(:,:),cf_ecosys%matrix_K(:,:)),cs_ecosys%Cpool(:)) 
                           
     cs_ecosys%Cpool(:) = cs_ecosys%Cpool(:) + cf_ecosys%dC_dt(:) * dt
                           
  end if
  
end subroutine update_ecosysC

subroutine matrix_diagnostic(cf_ecosys,cs_ecosys,SimTime)
  use EcosystemCarbonStateType, only : ecosystem_carbonstate_type
  use EcosystemCarbonFluxType,  only : ecosystem_carbonflux_type
  use SimulationTimeType, only : simulation_time_type
  use ParameterType, only : npools
  
  type(ecosystem_carbonflux_type) ,intent(inout) :: cf_ecosys
  type(ecosystem_carbonstate_type),intent(inout) :: cs_ecosys
  type(simulation_time_type)      ,intent(in)    :: SimTime
  
  
  real(kind=8),dimension(1:npools,1:npools)      :: matrix_AKinv
  real(kind=8),dimension(1:npools,1:npools)      :: matrix_AK
  real(kind=8) dt

  dt = SimTime%get_dtime()
  
  ! X(t) = -(AK)**(-1)BI(t) + (AK)**(-1) dX/dt 
  !      = Xc(t) - Xp(t)
  
  ! calculate "AK"
  matrix_AK = matmul(cf_ecosys%matrix_A(:,:),cf_ecosys%matrix_K(:,:))
  ! calculate "(AK)**-1"
  call inverse(matrix_AK,matrix_AKinv,npools)
  
  ! a) I = Cinput  
  cf_ecosys%diag_I = cf_ecosys%Cinput
  ! b) ResidenceTime
  ! ResidenceTime = AK**(-1)B
  cf_ecosys%diag_ResidenceTime = -matmul(matrix_AKinv,cf_ecosys%matrix_B)   ! Delete to let trainee complete
  
  ! c) C storage capacity
  ! Xc(t) = -(AK)**(-1)BI(t) * dt
  cf_ecosys%diag_Capacity      = -matmul(matrix_AKinv,cf_ecosys%matrix_B) * cf_ecosys%Cinput ! Delete to let trainee complete
  ! d) C storage potential
  ! Xp(t) = -(AK)**(-1) dX/dt * dt
  cf_ecosys%diag_Potential     = -matmul(matrix_AKinv,cf_ecosys%dC_dt) ! Delete to let trainee complete
  
end subroutine matrix_diagnostic
    
subroutine matrix_spinup(cf_ecosys,cs_ecosys,SimTime)
  use EcosystemCarbonStateType, only : ecosystem_carbonstate_type,&
          ifoliage,iwoody,imetlit,istrlit,imicsoil,islowsoil,ipasssoil
  use EcosystemCarbonFluxType,  only : ecosystem_carbonflux_type
  use SimulationTimeType
  use ParameterType, only : npools

  type(ecosystem_carbonflux_type) ,intent(inout) :: cf_ecosys
  type(ecosystem_carbonstate_type),intent(inout) :: cs_ecosys
  type(simulation_time_type)      ,intent(in)    :: SimTime
  
  real(kind=8),save,dimension(1:npools,1:npools) :: matrix_AK=0
  real(kind=8),dimension(1:npools,1:npools)      :: matrix_AKinv
  real(kind=8),save,dimension(1:npools)          :: matrix_BI=0
  integer,save :: nacc=0
  
  matrix_AK  = matrix_AK + matmul(cf_ecosys%matrix_A(:,:),cf_ecosys%matrix_K(:,:))
  matrix_BI  = matrix_BI + cf_ecosys%Cinput * cf_ecosys%matrix_B
              
  nacc = nacc + 1
  
  if(SimTime%is_endof_forcing())then
      
      matrix_AK = matrix_AK/nacc
      matrix_BI = matrix_BI/nacc
      
      call inverse(matrix_AK,matrix_AKinv,npools)
      cs_ecosys%Cpool = -matmul(matrix_AKinv,matrix_BI)

      nacc      = 0
      matrix_AK = 0
      matrix_BI = 0
      
   end if

end subroutine matrix_spinup
 
subroutine open_io_file(lu_Cinput,lu_param,lu_init,lu_out,file_Cinput,file_param,file_init,file_out)

  integer,intent(in) :: lu_Cinput
  integer,intent(in) :: lu_param
  integer,intent(in) :: lu_init
  integer,intent(in) :: lu_out

  character(len=150),intent(in) :: file_Cinput
  character(len=150),intent(in) :: file_param
  character(len=150),intent(in) :: file_init
  character(len=150),intent(in) :: file_out
  
  open(unit=lu_Cinput ,file=trim(adjustl(file_Cinput)))
  open(unit=lu_param  ,file=trim(adjustl(file_param)))
  open(unit=lu_init   ,file=trim(adjustl(file_init)))
  open(unit=lu_out    ,file=trim(adjustl(file_out)))
  write(lu_out,"(2A20,16A20)")'Year','Doy','Cinput(gC/s)','Foliage(gC/m2)','Woody(gC/m2)','Met_litter(gC/m2)','Str_litter(gC/m2)',&
                              'Fast_soil(gC/m2)','Slow_soil(gC/m2)','Pass_soil(gC/m2)','I(gC/m2/s)','Residence time(yr)',&
                              'Capacity(gC/m2)','Potential(gC/m2)'

end subroutine open_io_file

subroutine C_input(lu_Cinput,cf_ecosys)
  use EcosystemCarbonFluxType
  integer,intent(in) :: lu_Cinput
  type(ecosystem_carbonflux_type),intent(inout) :: cf_ecosys

  read(lu_Cinput,*)cf_ecosys%Cinput

end subroutine C_input

subroutine write_output(lu_out,cf_ecosys,cs_ecosys,SimTime,is_write_yearly,is_diagnostic)
  use EcosystemCarbonStateType
  use EcosystemCarbonFluxType
  use SimulationTimeType
  use ParameterType, only : npools

  integer,intent(in) :: lu_out
  type(ecosystem_carbonstate_type),intent(inout) :: cs_ecosys
  type(ecosystem_carbonflux_type) ,intent(inout) :: cf_ecosys
  type(simulation_time_type)      ,intent(in)    :: SimTime
  logical                         ,intent(in)    :: is_write_yearly
  logical                         ,intent(in)    :: is_diagnostic

  integer i
  integer(kind=8) year,doy
  logical is_write
  integer,save :: nave=0
  
  call SimTime%get_curr_time(year,doy)
  
  if(SimTime%is_startof_year())then
     nave = 0
     cf_ecosys%Cinput_out = 0.
     cs_ecosys%Cpool_out  = 0.
     if(is_diagnostic)then
        cf_ecosys%diag_I_out = 0.
        cf_ecosys%diag_ResidenceTime_out = 0.
        cf_ecosys%diag_Capacity_out = 0.
        cf_ecosys%diag_Potential_out = 0.
     end if
  end if
 
  if(is_write_yearly)then ! write output yearly
      if(SimTime%is_endof_year())then
          is_write = .True.
      else
          is_write = .False.
      end if
  else ! write output daily
      is_write = .True.
  end if
  
  nave = nave + 1
  cf_ecosys%Cinput_out   = cf_ecosys%Cinput_out   + cf_ecosys%Cinput
  cs_ecosys%Cpool_out(:) = cs_ecosys%Cpool_out(:) + cs_ecosys%Cpool(:)
  
  if(is_diagnostic)then
     cf_ecosys%diag_I_out             = cf_ecosys%diag_I_out             + cf_ecosys%diag_I
     cf_ecosys%diag_ResidenceTime_out = cf_ecosys%diag_ResidenceTime_out + sum(cf_ecosys%diag_ResidenceTime(:))  
     cf_ecosys%diag_Capacity_out      = cf_ecosys%diag_Capacity_out      + sum(cf_ecosys%diag_Capacity(:))  
     cf_ecosys%diag_Potential_out     = cf_ecosys%diag_Potential_out     + sum(cf_ecosys%diag_Potential(:))  
  end if
  
  if(is_write)then
     cf_ecosys%Cinput_out = cf_ecosys%Cinput_out / nave
     cs_ecosys%Cpool_out  = cs_ecosys%Cpool_out  / nave
     if(is_diagnostic)then
        cf_ecosys%diag_I_out             = cf_ecosys%diag_I_out             / nave * 3600 * 24 * 365
        cf_ecosys%diag_ResidenceTime_out = cf_ecosys%diag_ResidenceTime_out / nave / 3600 / 24 / 365
        cf_ecosys%diag_Capacity_out      = cf_ecosys%diag_Capacity_out      / nave 
        cf_ecosys%diag_Potential_out     = cf_ecosys%diag_Potential_out     / nave 
     end if
     
     write(lu_out,"(2I20,12E20.9)")year,doy,cf_ecosys%Cinput_out,(cs_ecosys%Cpool_out(i),i=1,npools),&
           cf_ecosys%diag_I_out,cf_ecosys%diag_ResidenceTime_out,cf_ecosys%diag_Capacity_out,cf_ecosys%diag_Potential_out



  end if
  
end subroutine write_output
  
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer,intent(in) :: n
real(kind=8),intent(in)  :: a(n,n)
real(kind=8),intent(out) :: c(n,n)
real(kind=8) :: L(n,n), U(n,n), aa(n,n), b(n), d(n), x(n)
real(kind=8) :: coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

aa=a

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=aa(i,k)/aa(k,k)
      L(i,k) = coeff
      do j=k+1,n
         aa(i,j) = aa(i,j)-coeff*aa(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = aa(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

    


