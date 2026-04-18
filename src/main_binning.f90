! Binning analysis for part 2
program main_binning
  use binning
  implicit none
  integer, parameter :: L = 20
  integer, parameter :: N = L * L
  integer(kind=8), parameter:: DISCARD_MCS = 1000_8
  character(len=*), parameter:: TS_FILE      = '../results/part2/timeseries_part2.dat'
  character(len=*), parameter:: SUMMARY_FILE = '../results/part2/summary_part2.dat'
  character(len=*), parameter:: BIN_FILE     = '../results/part2/binning_part2.dat'
  character(len=256):: line
  integer:: ios, unit_in, unit_sum, unit_bin
  integer:: k, maxk, q, nmax, nblocks, block_size, best_block
  integer(kind=8):: mcs
  integer, allocatable:: counts(:), fill(:)
  double precision:: temp, E
  double precision, allocatable:: temps(:), Eseries(:,:)
  double precision:: avg, sig
  double precision:: mean_e_spin, err_e_spin
  double precision:: mean_e, err_e
  logical:: first_block

  unit_in = 20
  unit_sum = 30
  unit_bin = 31

  ! First pass: find max temperature index
  maxk = 0
  open(unit_in, file=TS_FILE, status='old')
  read(unit_in, '(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, q, E
    if (ios /= 0) exit
    if (mcs >= DISCARD_MCS) then
      if (k > maxk) maxk = k
    endif
  enddo
  close(unit_in)

  if (maxk <= 0) stop 'No measurements found after thermalization cut.'

  allocate(counts(maxk), fill(maxk), temps(maxk))
  counts = 0
  fill = 0
  temps = 0.d0

  ! Second pass: count samples per temperature
  open(unit_in, file=TS_FILE, status='old')
  read(unit_in, '(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, q, E
    if (ios /= 0) exit
    if (mcs >= DISCARD_MCS) then
      counts(k) = counts(k) + 1
      temps(k) = temp
    endif
  enddo
  close(unit_in)

  nmax = maxval(counts)
  allocate(Eseries(nmax, maxk))
  Eseries = 0.d0

  ! Third pass: store energy per spin series
  open(unit_in, file=TS_FILE, status='old')
  read(unit_in, '(A)', iostat=ios) line
  do
    read(unit_in, *, iostat=ios) mcs, k, temp, q, E
    if (ios /= 0) exit
    if (mcs >= DISCARD_MCS) then
      fill(k) = fill(k) + 1
      Eseries(fill(k), k) = E / dble(N)
      temps(k) = temp
    endif
  enddo
  close(unit_in)

  ! Output files
  open(unit_sum, file=SUMMARY_FILE, status='replace')
  open(unit_bin, file=BIN_FILE, status='replace')

  write(unit_sum,'(A)') '# k T ndata meanE meanE_per_spin errE errE_per_spin best_block'
  write(unit_bin,'(A)') '# k T block_size nblocks meanE meanE_per_spin errE errE_per_spin'

  do k = 1, maxk
    if (counts(k) <= 1) cycle

    first_block = .true.
    best_block = 1
    mean_e_spin = 0.d0
    err_e_spin = 0.d0

    block_size = 1
    do while (counts(k) / block_size >= 2)
      nblocks = counts(k) / block_size

      call average_block(Eseries(1:counts(k), k), block_size, avg, sig)

      write(unit_bin,'(I4,1X,F12.6,1X,I8,1X,I8,1X,F20.10,1X,F20.10,1X,F20.10,1X,F20.10)') &
           k, temps(k), block_size, nblocks, avg*dble(N), avg, sig*dble(N), sig

      if (first_block) then
        mean_e_spin = avg
        err_e_spin = sig
        best_block = block_size
        first_block = .false.
      else
        ! Keep the largest error found
        if (sig > err_e_spin) then
          err_e_spin = sig
          best_block = block_size
        endif
      endif

      block_size = 2 * block_size
    enddo

    mean_e = mean_e_spin * dble(N)
    err_e  = err_e_spin  * dble(N)

    write(unit_sum,'(I4,1X,F12.6,1X,I8,1X,F20.10,1X,F20.10,1X,F20.10,1X,F20.10,1X,I8)') &
         k, temps(k), counts(k), mean_e, mean_e_spin, err_e, err_e_spin, best_block
  enddo

  close(unit_sum)
  close(unit_bin)

  deallocate(counts, fill, temps, Eseries)

end program main_binning