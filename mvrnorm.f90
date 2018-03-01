module mvrnorm

implicit none

! When this module is evoked these become global variables, at least under scope
!   of the main program. 
real (kind = 4) :: foo
real (kind = 8) :: bar
integer (kind = 4) :: n_stable_yrs, n_sims, init_stable_yr, end_stable_yr, hunt_yr0, n_yr_proj
integer (kind = 4) :: first_stable_yr, last_stable_yr, stable_len, n_lambda_yrs
integer ( kind = 4 ) :: abun_gap  ! number of yrs between last abundance estimate and init hunt yr
integer ( kind = 4 ) :: proj_len  ! number of yrs in projection (from last abundance estimate to end of waiver)
real (kind = 8) :: m_hunt_mu, ci_1, ci_2
logical :: use_cov, print_inits, first_poisson
integer ( kind = 4 ) :: seed

real ( kind = 8 ), allocatable, dimension ( :, : ) :: varcov         ! variance covariance matrix passed as input to mvrnorm
real ( kind = 8 ), allocatable, dimension ( : ) :: mu                ! vector of expected value passed as input to mvrnorm 
real ( kind = 8 ), allocatable, dimension ( :, : ) :: rand_stable    ! vector of deviates returned from mvrnorm -- bootstrap PCFG abundance during stable period
real ( kind = 8 ), allocatable, dimension ( : ) :: mu_hunt           ! expected average annual PCFG whaling mortality (1.6 per yr)
real ( kind = 8 ), allocatable, dimension ( : , : ) :: N_proj_pcfg   ! matrix to hold future projections
real ( kind = 8 ), allocatable, dimension ( : ) :: lambda_boot       ! annual rates of change from the bootstrap sample
real ( kind = 8 ), allocatable, dimension ( : ) :: mu_lambda_boot    ! mean log(lambda) for the boostrap sample
real ( kind = 8 ), allocatable, dimension ( : ) :: sd_lambda_boot    ! sd for log(lambda) for the bootstrap
real ( kind = 8 ), allocatable, dimension ( : ) :: mu_r              ! mean r for stable period: 2002-2015
real ( kind = 8 ), allocatable, dimension ( : ) :: lambda_hunt_proj  ! annual lambdas for the projected population sizes after the hunt
real ( kind = 8 ), allocatable, dimension ( : ) :: log_lambda_hunt_proj  ! log of lambdas for the projected population sizes after the hunt
real ( kind = 8 ), allocatable, dimension ( : ) :: hunt_m     ! vector of realized hunting mortality by year
real ( kind = 8 ), allocatable, dimension ( : ) :: normal_rv  ! vector of normal random variates for lambdas by year (changes in N assumed log-normal)
real ( kind = 8 ), allocatable, dimension ( : ) :: mu_r_proj  ! mean r (across years) for the PROJECTED population sizes after the hunt (for boostrap dataset i)
! real ( kind = 8 ), allocatable, dimension ( : ) :: N_rand  ! vector of random deviates from multivariate normal abundance estimates 
! Not sure these variables are needed here 
! REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
!                       vsmall = TINY(1.0), vlarge = HUGE(1.0)
PRIVATE            :: integral
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

CONTAINS

!****************************************************************************JRB
function mean(x) result(fn_val)
  ! Calculate arithmetic mean of vector of reals
  real (kind = 8) :: fn_val
  real (kind = 8) :: x(:)
  fn_val = sum(x) / size(x)
  return
end function mean
!****************************************************************************JRB
function var(x) result(fn_val)
  ! Calculate variance of vector of reals
  real (kind = 8) :: mu, fn_val
  real (kind = 8) :: x(:)
  mu = mean(x)
  fn_val = sum((mu - x)**2) / (size(x) - 1)
  return
end function var  
!****************************************************************************JRB
function sd(x) result(fn_val)
  ! Calculate standard deviation of vector of reals
  real (kind = 8) :: fn_val
  real (kind = 8) :: x(:)
  fn_val = sqrt(var(x))
  return
end function sd
!****************************************************************************JRB
! TODO: Move this subroutine into a separate module -- `noaa_sla.f90`
subroutine read_inits()
  !  Initialize variables for NOAA advection-diffusion SLA 
  open(unit = 13, file = 'config.ini')
  read(13, *) n_stable_yrs    ! 14     # number of years during 'stable' period (i.e. dimensions of covariance matrix): 2002-2015
  read(13, *) first_stable_yr ! 2002
  read(13, *) last_stable_yr  ! 2015
  read(13, *) n_sims          ! 50000  # number of random variates to generate, ie number of bootstrap samples
  ! print *, 'n_sims', n_sims
  ! print *, "STOPPING"
  ! stop
  read(13, *) seed            ! 123    # random number seed
  read(13, *) init_stable_yr  ! 7      # index of first year (2002) in 'stable' period
  read(13, *) end_stable_yr   ! 20     # index of last year (2015) in 'stable' period
  read(13, *) hunt_yr0        ! 2019   # when will the hunt begin (hypothetical)?
  read(13, *) n_yr_proj       ! 10     # number of years to project population forward starting with the first HUNT year
  read(13, *) use_cov         ! .true. # use variance covariance matrix?  boolean
  read(13, *) m_hunt_mu       ! 1.6    # mean annual hunting mortality once the hunt begins
  read(13, *) ci_1  ! 0.60   # confidence limits for projections. e.g 0.9 = 90% CI, so lower limit is 5th percentile.  Use 0.60 to get the 20th percentile.
  read(13, *) ci_2  ! 0.80   # This should a the wider level (e.g., higher clim2 > clim1)
  read(13, *) print_inits     ! print out these initialized variables for checking?
  close(13)
  if (print_inits) then
    print *, 'inits have been read'
    print *, 'first_stable_yr', first_stable_yr
    print *, 'last_stable_yr', last_stable_yr  
    print *, 'seed', seed
    print *, 'init_stable_yr', init_stable_yr
    print *, 'end_stable_yr', end_stable_yr
    print *, 'hunt_yr0', hunt_yr0
    print *, 'n_yr_proj', n_yr_proj
    print *, 'n_sims', n_sims
    print *, 'use_cov', use_cov
    print *, 'm_hunt_mu', m_hunt_mu
    print *, 'ci_1', ci_1  
    print *, 'ci_2', ci_2
  end if
  return  
end subroutine read_inits

subroutine allocate_arrays()
!*****************************************************************************80
  integer :: i
! initialize time windows for simulations
  abun_gap = hunt_yr0 - last_stable_yr - 1
  proj_len = abun_gap + n_yr_proj
  stable_len = last_stable_yr - first_stable_yr + 1
  n_lambda_yrs = stable_len - 1
! Initialize arrays for bootstrap and forward projections
  allocate ( lambda_boot(1:(stable_len - 1)) )  ! annual rates of change from the bootstrap sample
  ! allocate ( N_rand(stable_len) )       ! vector of random deviates from multivariate normal abundance estimates during stable period
  allocate ( mu_lambda_boot(1:n_sims) )   ! mean of log(lambda_boot) for a bootstrap sample
  allocate ( sd_lambda_boot(1:n_sims) )   ! SD of log(lambda_boot) for a bootstrap sample
  allocate ( N_proj_pcfg(1:n_sims, 1:proj_len) )  ! matrix to hold future projections
  ! print *, 'n_sims', n_sims
  ! print *, 'shape(N_proj_pcfg)', shape(N_proj_pcfg)
  ! print *, "STOPPING"
  ! stop      
  allocate ( mu_r(1:n_sims) )
  allocate ( mu_r_proj(1:n_sims) )
  allocate ( lambda_hunt_proj(1:proj_len) )
  allocate ( log_lambda_hunt_proj(1:proj_len) )
  allocate ( hunt_m(1:proj_len) )
  allocate ( normal_rv(1:proj_len) )
  allocate ( mu_hunt(1:proj_len) )
  allocate ( mu(1:n_stable_yrs) )
! Assign expectation for future PCFG hunting mortality 
  do i = 1, abun_gap
    mu_hunt(i) = 0.d0
  end do
  do i = (abun_gap + 1), proj_len
    mu_hunt(i) = m_hunt_mu
  end do
!  Read point estimates of abundance during 'stable period'
  open(unit = 12, file = 'mu.dat')
  read(12, *) mu
  close(12)
! Read variance-covariance matrix of abundance estimates during 'stable period'
  allocate ( varcov(1:n_stable_yrs,1:n_stable_yrs) )
  open(unit = 11, file = 'cov.dat')
  read(11, *) varcov
  close(11)
return
end subroutine allocate_arrays

subroutine get_unit ( iunit )
!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) :: i
  integer ( kind = 4 ) :: ios
  integer ( kind = 4 ) :: iunit
  logical :: lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end subroutine get_unit

subroutine multinormal_sample ( m, n, a, mu, seed, x )
!*****************************************************************************80
!
!! MULTINORMAL_SAMPLE samples a multivariate normal distribution.
!
!  Discussion:
!
!    The multivariate normal distribution for the M dimensional vector X
!    has the form:
!
!      pdf(X) = (2*pi*det(A))**(-M/2) * exp(-0.5*(X-MU)'*inverse(A)*(X-MU))
!
!    where MU is the mean vector, and A is a positive definite symmetric
!    matrix called the variance-covariance matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the space.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) A(M,M), the variance-covariance 
!    matrix.  A must be positive definite symmetric.
!
!    Input, real ( kind = 8 ) MU(M), the mean vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Output, real ( kind = 8 ) X(M), the points.
!
  implicit none

  integer ( kind = 4 ) :: m
  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: a(m,m)
  integer ( kind = 4 ) :: info
  integer ( kind = 4 ) :: j
  real ( kind = 8 ) :: mu(m)
  real ( kind = 8 ) :: r(m,m)
  integer ( kind = 4 ) :: seed
  real ( kind = 8 ) :: x(m,n)
!
!  Compute the upper triangular Cholesky factor R of the variance-covariance
!  matrix.
!
  r(1:m,1:m) = a(1:m,1:m)

  call r8po_fa ( m, r, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MULTINORMAL_SAMPLE - Fatal error!'
    write ( *, '(a)' ) &
      '  The variance-covariance matrix is not positive definite symmetric.'
    stop
  end if
!
!  Get an MxN matrix of samples of the 1D normal distribution with mean 0
!  and variance 1.  
!
  call r8vec_normal_01 ( m*n, seed, x(1:m,1:n) )
!
!  Compute R' * X.
!  We actually carry out this computation in the equivalent form X' * R.
!
  do j = 1, n
    x(1:m,j) = mu(1:m) + matmul ( x(1:m,j), r(1:m,1:m) )
  end do

  return
end subroutine multinormal_sample

function r8_uniform_01 ( seed )
!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) :: k
  real ( kind = 8 ) :: r8_uniform_01
  integer ( kind = 4 ) :: seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01

subroutine r8mat_print ( m, n, a, title )
!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) :: m
  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: a(m,n)
  character ( len = * ) :: title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end subroutine r8mat_print

subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )
!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) :: m
  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: a(m,n)
  character ( len = 14 ) :: ctemp(incx)
  integer ( kind = 4 ) :: i
  integer ( kind = 4 ) :: i2hi
  integer ( kind = 4 ) :: i2lo
  integer ( kind = 4 ) :: ihi
  integer ( kind = 4 ) :: ilo
  integer ( kind = 4 ) :: inc
  integer ( kind = 4 ) :: j
  integer ( kind = 4 ) :: j2
  integer ( kind = 4 ) :: j2hi
  integer ( kind = 4 ) :: j2lo
  integer ( kind = 4 ) :: jhi
  integer ( kind = 4 ) :: jlo
  character ( len = * ) :: title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end subroutine r8mat_print_some

subroutine r8mat_write ( output_filename, m, n, table )
!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) :: m
  integer ( kind = 4 ) :: n

  integer ( kind = 4 ) :: j
  character ( len = * ) :: output_filename
  integer ( kind = 4 ) :: output_status
  integer ( kind = 4 ) :: output_unit
  character ( len = 30 ) :: string
  real ( kind = 8 ) :: table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end subroutine r8mat_write

subroutine r8po_fa ( n, a, info )
!*****************************************************************************80
!
!! R8PO_FA factors an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!    The positive definite symmetric matrix A has a Cholesky factorization
!    of the form:
!
!      A = R' * R
!
!    where R is an upper triangular matrix with positive elements on
!    its diagonal.  This routine overwrites the matrix A with its
!    factor R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the matrix in R8PO storage.
!    On output, the Cholesky factor R in R8GE storage.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal return.
!    K, error condition.  The principal minor of order K is not
!    positive definite, and the factorization was not completed.
!
  implicit none

  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: a(n,n)
  integer ( kind = 4 ) :: i
  integer ( kind = 4 ) :: info
  integer ( kind = 4 ) :: j
  integer ( kind = 4 ) :: k
  real ( kind = 8 ) :: s

  do j = 1, n

    do k = 1, j - 1
      a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
    end do

    s = a(j,j) - sum ( a(1:j-1,j)**2 )

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( s )

  end do

  info = 0
!
!  Since the Cholesky factor is stored in R8GE format, be sure to
!  zero out the lower triangle.
!
  do i = 1, n
    do j = 1, i-1
      a(i,j) = 0.0D+00
    end do
  end do

  return
end subroutine r8po_fa

SUBROUTINE integral(a, b, result, dk)
!*****************************************************************************80  
!     Gaussian integration of exp(k.cosx) from a to b.
REAL (dp), INTENT(IN) :: dk
REAL, INTENT(IN)      :: a, b
REAL, INTENT(OUT)     :: result

!     Local variables

REAL (dp)  :: xmid, range, x1, x2,                                    &
  x(3) = (/0.238619186083197_dp, 0.661209386466265_dp, 0.932469514203152_dp/), &
  w(3) = (/0.467913934572691_dp, 0.360761573048139_dp, 0.171324492379170_dp/)
INTEGER    :: i

xmid = (a + b)/2._dp
range = (b - a)/2._dp

result = 0._dp
DO i = 1, 3
  x1 = xmid + x(i)*range
  x2 = xmid - x(i)*range
  result = result + w(i)*(EXP(dk*COS(x1)) + EXP(dk*COS(x2)))
END DO

result = result * range
RETURN
END SUBROUTINE integral

subroutine r8vec_normal_01 ( n, seed, x )
!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is 
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) :: n

  integer ( kind = 4 ) :: m
  integer ( kind = 4 ), save :: made = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r(n+1)
!  real ( kind = 8 ) :: r8_uniform_01  ! declaration here (r8_uniform_01 is a function) was confusing compiler
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) :: seed
  real ( kind = 8 ) :: x(n)
  integer ( kind = 4 ) :: x_hi_index
  integer ( kind = 4 ) :: x_lo_index
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2 * m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2 * m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end subroutine r8vec_normal_01

subroutine r8vec_print ( n, a, title )
!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: a(n)
  integer ( kind = 4 ) :: i
  character ( len = * ) :: title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end subroutine r8vec_print

subroutine r8vec_uniform_01 ( n, seed, r )
!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) :: n

  integer ( kind = 4 ) :: i
  integer ( kind = 4 ) :: k
  integer ( kind = 4 ) :: seed
  real ( kind = 8 ) :: r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end subroutine r8vec_uniform_01

subroutine timestamp ( )
!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) :: ampm
  integer ( kind = 4 ) :: d
  integer ( kind = 4 ) :: h
  integer ( kind = 4 ) :: m
  integer ( kind = 4 ) :: mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp

end module mvrnorm