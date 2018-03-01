program main
!  Author: John R. Brandon
!  Purpose: Implement NOAA advection-diffusion gray whale SLA for Makah whaling  

  use mvrnorm  ! For generating multivariate random normal deviates
  use random   ! Various probability distributions
!
  implicit none

  integer ( kind = 4 ) :: i, j
  character ( len = 255 ) :: output_filename
  character ( len = 255 ) :: string

  call read_inits()  ! read initial values (routine is in module `mvrnorm.f90`)
  call allocate_arrays()  ! initial values include dimensions and some values for arrays

!  Generate the random deviates: NOTE each row is a bootstrap sample, with columns corresponding to years
  allocate ( rand_stable(1:n_stable_yrs, 1:n_sims) )  ! m = n_stable_yrs; n = n_sims
  call multinormal_sample ( n_stable_yrs, n_sims, varcov, mu, seed, rand_stable )
  
  ! write (*,*) 'n_stable_yrs', n_stable_yrs  
  ! print *, 'shape(rand_stable)', shape(rand_stable)
  ! print *, 'proj_len', proj_len
  ! do j = 1, n_stable_yrs
  !   write(*,*) j, 'rand_stable(j, 2)', rand_stable(j, 2)
  ! end do

  ! first_poisson = .false.

  do j = 1, n_sims

    ! Find annual rates of change from the bootstrap sample
    do i = 2, n_stable_yrs
      lambda_boot(i - 1) = rand_stable(i, j) / rand_stable(i - 1, j)
    end do

    ! do i = 1, (n_stable_yrs - 1)
    !   print *, 'lambda_boot(i)', lambda_boot(i)
    ! end do 

    ! Find mean and SD of log(lambda) for this boostrap sample
    mu_lambda_boot(j) = mean(log(lambda_boot))
    sd_lambda_boot(j) = sd(log(lambda_boot))
    ! print *, 'mu_lambda_boot', mu_lambda_boot(1), 'sd_lambda_boot', sd_lambda_boot(1)

    ! Project N for first future year
    ! Take last abundance estimate and generate log-normal change in abundance conditional on lambdas  
    normal_rv(1) = random_normal(mu_lambda_boot(j), sd_lambda_boot(j))    ! generate normal random variate for abundance | lambda
    hunt_m(1) = random_Poisson(mu_hunt(1), .true.)
    N_proj_pcfg(j, 1) = rand_stable(n_stable_yrs, j) * exp(normal_rv(1))  ! log-normal change in abundance     
    N_proj_pcfg(j, 1) = N_proj_pcfg(j, 1) - hunt_m(1)  ! subtract future expected hunt mortality
    ! N_proj_pcfg(j, 1) = N_proj_pcfg(j, 1) - random_Poisson(mu_hunt(1), .true.)  ! subtract future expected hunt mortality
    ! N_proj_pcfg(j, 1) = N_proj_pcfg(j, 1) - mu_hunt(1)  ! subtract future expected hunt mortality    
    
    ! Record lambda for first year of projection
    lambda_hunt_proj(1) = N_proj_pcfg(j, 1) / rand_stable(n_stable_yrs, j)

    ! # Project N for subsequent future years  
    do i = 2, proj_len
      normal_rv(i) = random_normal(mu_lambda_boot(j), sd_lambda_boot(j))  ! generate normal random variate for abundance | lambda
      hunt_m(i) = random_Poisson(mu_hunt(i), .true.)      
      N_proj_pcfg(j, i) = N_proj_pcfg(j, i - 1) * exp(normal_rv(i))       ! log-normal change in abundance     
      N_proj_pcfg(j, i) = N_proj_pcfg(j, i) - hunt_m(i)  ! subtract future expected hunt mortality
      ! N_proj_pcfg(j, i) = N_proj_pcfg(j, i) - random_Poisson(mu_hunt(i), .true.)  ! subtract future expected hunt mortality
      ! N_proj_pcfg(j, i) = N_proj_pcfg(j, i) - mu_hunt(i)  ! subtract future expected hunt mortality
      lambda_hunt_proj(i) = N_proj_pcfg(j, i) / N_proj_pcfg(j, i - 1)  ! Find the yearly lambdas for the projected abundance estimates
    end do
  
    ! find the mean r (across years) for the PROJECTED population sizes after the hunt (for boostrap dataset i)
    log_lambda_hunt_proj = log(lambda_hunt_proj) 
    mu_r_proj(j) = mean(log_lambda_hunt_proj((hunt_yr0 - last_stable_yr):proj_len))

    ! do i = 1, proj_len  ! (hunt_yr0 - last_stable_yr), proj_len
    !   print *, i, normal_rv(i), N_proj_pcfg(j, i), lambda_hunt_proj(i), log_lambda_hunt_proj(i), hunt_m(i), mu_hunt(i) ! , random_Poisson(mu_hunt(i), .true.)
    ! end do
    
    ! print *, j, 'mu_r_proj(j)', mu_r_proj(j)
    
    ! write(*,*) (N_proj_pcfg(j, i), lambda_hunt_proj(i), i = 1, proj_len)
    ! print *, log(lambda_hunt_proj((hunt_yr0 - last_stable_yr):proj_len))
    
    ! mu_r_proj(proj_len) = 0.d0
    ! do i = 1, proj_len
    !   write (*, *) i, 'mu_r_proj(i)', mu_r_proj(i)
    ! end do
     
    ! # this is the mean r for 2002-2015
    ! mu.r[i] = mu.mean[i]

  end do         

  ! call r8mat_print ( n_stable_yrs, n_sims, rand_stable, '  Multinormal variates:' )
  ! call r8mat_print ( n_sims, proj_len, N_proj_pcfg, '  N_proj_pcfg:' )

!  Write the data to a file.
  ! write ( output_filename, '(a, i2.2, a, i5.5, a)' ) 'normal_', m, '_', n, '.txt'
  output_filename = 'mvnorm.dat'
  call r8mat_write ( output_filename, n_stable_yrs, n_sims, rand_stable )

  output_filename = 'N_proj.dat'
  ! call r8mat_write ( output_filename, n_sims, proj_len, N_proj_pcfg )  
  call r8mat_write ( output_filename, proj_len, n_sims, transpose(N_proj_pcfg) )  
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to the file "' &
     // trim ( output_filename ) // '".'

  !  Terminate.
  ! call timestamp ( )

  stop
end program main

