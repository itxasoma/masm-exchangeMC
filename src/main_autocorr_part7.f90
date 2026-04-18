program main_autocorr_part7
  implicit none
  character(len=256):: infile, outfile, observable, arg
  character(len=256):: line
  integer:: nargs, ios
  integer:: nrows, i, j, lag, maxlag
  integer(kind=8), allocatable:: mcs_all(:)
  integer, allocatable:: k_all(:), q_all(:)
  double precision, allocatable:: t_all(:), e_all(:)
  double precision, allocatable:: x(:), rho(:), tau_run(:)
  double precision:: temp_target, temp_tol, xmean, var0, tau_int
  logical:: select_temp
  integer:: nsel
  integer(kind=8):: mcs_tmp
  integer:: k_tmp, q_tmp
  double precision:: t_tmp, e_tmp

  ! Co-pilot suggestion for selecting and checking the arguments
  nargs = command_argument_count()
  if (nargs < 4) then
    print *, 'Usage: ./main_autocorr_part7.x input.dat output.dat observable[q|e] maxlag [temperature]'
    stop
  endif
  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)
  call get_command_argument(3, observable)
  call get_command_argument(4, arg)
  read(arg, *) maxlag

  ! Co-pilot suggestion for selecting the temperature
  select_temp = .false.
  temp_target = 0.d0
  temp_tol = 1.d-8
  if (nargs >= 5) then
    call get_command_argument(5, arg)
    read(arg, *) temp_target
    select_temp = .true.
  endif

  ! Read the input file to count the number of rows and allocate arrays
  open(10, file=trim(infile), status='old')
    nrows = 0
    read(10, '(A)', iostat=ios) line
    do
        read(10, *, iostat=ios) mcs_tmp, k_tmp, t_tmp, q_tmp, e_tmp
        if (ios /= 0) exit
        nrows = nrows + 1
    enddo
  close(10)
  allocate(mcs_all(nrows), k_all(nrows), q_all(nrows), t_all(nrows), e_all(nrows))

  ! Then, read the data:
  open(10, file=trim(infile), status='old')
    read(10, '(A)', iostat=ios) line
    do i = 1, nrows
        read(10, *, iostat=ios) mcs_all(i), k_all(i), t_all(i), q_all(i), e_all(i)
        if (ios /= 0) stop 'Error while reading input file.'
    enddo
  close(10)

  ! Choosing the data for the selected temperature
  if (select_temp) then
    nsel = 0
    do i = 1, nrows
      if (abs(t_all(i) - temp_target) < temp_tol) nsel = nsel + 1
    enddo
  else
    nsel = nrows
  endif
  allocate(x(nsel))

  ! Separate observable: q or e?
  j = 0
  if (trim(observable) == 'q') then
    do i = 1, nrows
      if ((.not. select_temp) .or. (abs(t_all(i) - temp_target) < temp_tol)) then
        j = j + 1
        x(j) = dble(q_all(i))
      endif
    enddo
  else if (trim(observable) == 'e') then
    do i = 1, nrows
      if ((.not. select_temp) .or. (abs(t_all(i) - temp_target) < temp_tol)) then
        j = j + 1
        x(j) = e_all(i)
      endif
    enddo
  else
    stop 'Observable must be q or e.'
  endif

  if (maxlag >= nsel) maxlag = nsel - 1
  if (maxlag < 1) stop 'maxlag must be at least 1.'

  ! Compute autocorrelation function and integrated autocorrelation time:
  allocate(rho(0:maxlag), tau_run(0:maxlag))

  xmean = sum(x) / dble(nsel)
  var0 = sum((x - xmean)**2) / dble(nsel)

  if (var0 <= 0.d0) stop 'Variance is zero.'

  do lag = 0, maxlag
    rho(lag) = 0.d0
    do i = 1, nsel - lag
      rho(lag) = rho(lag) + (x(i) - xmean) * (x(i+lag) - xmean)
    enddo
    rho(lag) = rho(lag) / dble(nsel - lag)
    rho(lag) = rho(lag) / var0
  enddo

  tau_int = 0.5d0
  tau_run(0) = tau_int

  do lag = 1, maxlag
    if (rho(lag) > 0.d0) then
      tau_int = tau_int + rho(lag)
    endif
    tau_run(lag) = tau_int
  enddo

  ! Write the results in output file
  open(20, file=trim(outfile), status='replace')
  write(20,'(A)') '# lag rho tau_int_running'
  do lag = 0, maxlag
    write(20,'(I8,1X,F18.10,1X,F18.10)') lag, rho(lag), tau_run(lag)
  enddo
  close(20)

  print *, 'Observable = ', trim(observable)
  if (select_temp) print *, 'Temperature = ', temp_target
  print *, 'Points used = ', nsel
  print *, 'Estimated tau_int = ', tau_int

  deallocate(mcs_all, k_all, q_all, t_all, e_all)
  deallocate(x, rho, tau_run)

end program main_autocorr_part7