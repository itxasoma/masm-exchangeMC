! BINNING FOR PART 4

program main_binning_part4
  use binning
  implicit none

  integer, parameter:: L = 8
  integer, parameter:: D = 3
  integer, parameter:: N = L**D    ! 512

  ! Time intervals: [lo, hi) MCS
  integer, parameter:: NINT = 4
  integer(kind=8), parameter:: INT_LO(NINT) = &
       [10000_8,  30000_8, 100000_8,  300000_8]
  integer(kind=8), parameter:: INT_HI(NINT) = &
       [30000_8, 100000_8, 300000_8, 1000001_8]

  character(len=*), parameter:: TS_FILE  = '../results/timeseries_part4.dat'
  character(len=*), parameter:: OUT_FILE = '../results/summary_part4.dat'

  integer:: ios, maxk, k, ii, block_size, nblocks
  integer(kind=8):: mcs
  integer:: Qi
  double precision:: temp, Eavg
  character(len=256):: line

  ! Per-temperature, per-interval counters and data
  integer, allocatable:: counts(:,:)   ! counts(k, interval)
  integer, allocatable:: fill(:,:)
  double precision, allocatable:: temps(:)
  double precision, allocatable:: Eseries(:,:,:)   ! (sample, k, interval)
  double precision, allocatable:: Q2series(:,:,:)  ! (sample, k, interval)

  integer:: nmax
  double precision:: avg_e, sig_e, avg_q2, sig_q2
  double precision:: best_sig_e, best_sig_q2
  integer:: unit_in, unit_out

  unit_in  = 20
  unit_out = 30

  ! --- Pass 1: find maxk and count samples per (k, interval) ---
  maxk = 0
  open(unit_in, file=TS_FILE, status='old')
  read(unit_in,'(A)', iostat=ios) line   ! skip header
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, Qi, Eavg
    if (ios /= 0) exit
    if (k > maxk) maxk = k
  enddo
  close(unit_in)

  allocate(counts(maxk, NINT), fill(maxk, NINT), temps(maxk))
  counts = 0;  fill = 0;  temps = 0.d0

  open(unit_in, file=TS_FILE, status='old')
  read(unit_in,'(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, Qi, Eavg
    if (ios /= 0) exit
    temps(k) = temp
    do ii = 1, NINT
      if (mcs >= INT_LO(ii) .and. mcs < INT_HI(ii)) then
        counts(k, ii) = counts(k, ii) + 1
      endif
    enddo
  enddo
  close(unit_in)

  nmax = maxval(counts)
  if (nmax < 2) stop 'Not enough data in time intervals'

  allocate(Eseries(nmax, maxk, NINT))
  allocate(Q2series(nmax, maxk, NINT))
  Eseries  = 0.d0
  Q2series = 0.d0

  ! --- Pass 2: store E/N and q^2 per interval ---
  open(unit_in, file=TS_FILE, status='old')
  read(unit_in,'(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, Qi, Eavg
    if (ios /= 0) exit
    do ii = 1, NINT
      if (mcs >= INT_LO(ii) .and. mcs < INT_HI(ii)) then
        fill(k, ii) = fill(k, ii) + 1
        Eseries(fill(k,ii),  k, ii) = Eavg / dble(N)
        Q2series(fill(k,ii), k, ii) = (dble(Qi) / dble(N))**2
      endif
    enddo
  enddo
  close(unit_in)

  ! --- Output ---
  open(unit_out, file=OUT_FILE, status='replace')
  write(unit_out,'(A)') &
    '# interval  k  T  ndata  meanE  errE  meanQ2  errQ2'

  do ii = 1, NINT
    do k = 1, maxk
      if (counts(k, ii) < 2) cycle

      ! Binning for E
      best_sig_e = 0.d0
      block_size = 1
      do while (counts(k,ii) / block_size >= 2)
        nblocks = counts(k,ii) / block_size
        call average_block(Eseries(1:counts(k,ii), k, ii), block_size, avg_e, sig_e)
        if (sig_e > best_sig_e) best_sig_e = sig_e
        block_size = 2 * block_size
      enddo

      ! Binning for Q2
      avg_e = sum(Eseries(1:counts(k,ii), k, ii)) / dble(counts(k,ii))
      avg_q2 = sum(Q2series(1:counts(k,ii), k, ii)) / dble(counts(k,ii))
      best_sig_q2 = 0.d0
      block_size = 1
      do while (counts(k,ii) / block_size >= 2)
        nblocks = counts(k,ii) / block_size
        call average_block(Q2series(1:counts(k,ii), k, ii), block_size, avg_q2, sig_q2)
        if (sig_q2 > best_sig_q2) best_sig_q2 = sig_q2
        block_size = 2 * block_size
      enddo

      write(unit_out,'(I2,1X,I4,1X,F10.5,1X,I8,4(1X,F16.8))') &
           ii, k, temps(k), counts(k,ii), avg_e, best_sig_e, avg_q2, best_sig_q2
    enddo
  enddo

  close(unit_out)
  print *, 'Written: ', OUT_FILE

  deallocate(counts, fill, temps, Eseries, Q2series)

end program main_binning_part4