! BINNING FOR PART 4

program main_binning_part4
  use binning
  implicit none

  integer, parameter:: L = 8
  integer, parameter:: D = 3
  integer, parameter:: N = L**D    ! 512

  ! Time intervals: [lo, hi) MCS
  integer, parameter:: NINT = 4
  integer(kind=8), parameter:: INT_LO(NINT) = [10000_8,  30000_8, 100000_8,  300000_8]
  integer(kind=8), parameter:: INT_HI(NINT) = [30000_8, 100000_8, 300000_8, 1000001_8]

  character(len=256):: ts_file, out_file, line
  integer:: ios, unit_in, unit_out
  integer:: k, maxk, q, nmax, block_size, interval
  integer(kind=8):: mcs
  integer, allocatable:: counts(:,:), fill(:,:)
  double precision:: temp, E
  double precision, allocatable:: temps(:)
  double precision, allocatable:: q2series(:,:,:), eseries(:,:,:)
  double precision:: avg_q2, err_q2, avg_e, err_e, best_err_q2, best_err_e
  integer:: nargs

  nargs = command_argument_count()
  if (nargs < 2) then
    print *, 'Usage: ./main_binning_part4.x input_timeseries.dat output_summary.dat'
    stop
  endif

  call get_command_argument(1, ts_file)
  call get_command_argument(2, out_file)

  unit_in = 20
  unit_out = 30

  ! Pass 1: find max temperature index
  maxk = 0
  open(unit_in, file=trim(ts_file), status='old')
  read(unit_in, '(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, q, E
    if (ios /= 0) exit
    if (k > maxk) maxk = k
  end do
  close(unit_in)

  if (maxk <= 0) stop 'No data found in input file.'

  allocate(counts(maxk, NINT), fill(maxk, NINT), temps(maxk))
  counts = 0
  fill = 0
  temps = 0.d0

  ! Pass 2: count samples in each interval for each temperature
  open(unit_in, file=trim(ts_file), status='old')
  read(unit_in, '(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, q, E
    if (ios /= 0) exit
    temps(k) = temp
    do interval = 1, NINT
      if (mcs >= INT_LO(interval) .and. mcs < INT_HI(interval)) then
        counts(k, interval) = counts(k, interval) + 1
      endif
    end do
  end do
  close(unit_in)

  nmax = maxval(counts)
  if (nmax <= 1) stop 'Not enough samples in selected intervals.'

  allocate(q2series(nmax, maxk, NINT))
  allocate(eseries(nmax, maxk, NINT))
  q2series = 0.d0
  eseries = 0.d0

  ! Pass 3: store q^2 and E/N
  open(unit_in, file=trim(ts_file), status='old')
  read(unit_in, '(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, q, E
    if (ios /= 0) exit
    do interval = 1, NINT
      if (mcs >= INT_LO(interval) .and. mcs < INT_HI(interval)) then
        fill(k, interval) = fill(k, interval) + 1
        q2series(fill(k, interval), k, interval) = (dble(q) / dble(N))**2
        eseries(fill(k, interval), k, interval)  = E / dble(N)
      endif
    end do
  end do
  close(unit_in)

  open(unit_out, file=trim(out_file), status='replace')
  write(unit_out,'(A)') '# interval k T ndata mean_q2 err_q2 mean_e_per_spin err_e_per_spin'

  do interval = 1, NINT
    do k = 1, maxk
      if (counts(k, interval) <= 1) cycle

      avg_q2 = sum(q2series(1:counts(k, interval), k, interval)) / dble(counts(k, interval))
      avg_e  = sum(eseries(1:counts(k, interval), k, interval)) / dble(counts(k, interval))

      best_err_q2 = 0.d0
      best_err_e  = 0.d0

      block_size = 1
      do while (counts(k, interval) / block_size >= 2)
        call average_block(q2series(1:counts(k, interval), k, interval), block_size, avg_q2, err_q2)
        call average_block(eseries(1:counts(k, interval), k, interval), block_size, avg_e, err_e)

        if (err_q2 > best_err_q2) best_err_q2 = err_q2
        if (err_e  > best_err_e ) best_err_e  = err_e

        block_size = 2 * block_size
      end do

      write(unit_out,'(I3,1X,I4,1X,F12.6,1X,I8,1X,F18.10,1X,F18.10,1X,F18.10,1X,F18.10)') &
           interval, k, temps(k), counts(k, interval), avg_q2, best_err_q2, avg_e, best_err_e
    end do
  end do

  close(unit_out)
  deallocate(counts, fill, temps, q2series, eseries)

end program main_binning_part4