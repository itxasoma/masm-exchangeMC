program main
  use rng_wrapper
  use GLOBAL
  use lattice
  use bonds
  use exchangeMC
  implicit none

  integer:: r, parity
  integer, allocatable:: spinsA(:,:), spinsB(:,:)
  integer, allocatable:: repA(:), repB(:), tempA(:), tempB(:)
  integer(kind=8):: mcs
  integer(kind=8), allocatable:: accA(:), attA(:), accB(:), attB(:)
  double precision, allocatable:: EA(:), EB(:)
  character(len=200):: inputfile, tsfile, swfile, tag
  ! Check post-transient swaps
  integer:: k
  integer(kind=8):: mcs_therm
  integer(kind=8), allocatable:: accA_eq(:), attA_eq(:), accB_eq(:), attB_eq(:)
  integer(kind=8), allocatable:: accA_old(:), attA_old(:), accB_old(:), attB_old(:)
  double precision:: rateA, rateB
  logical:: ok_all

  ! HISTOGRAM (on the run)
  integer, parameter:: NBINS = 2*64 + 1 ! 2N+1 bins for q = Q/N in [-1,1]
  integer(kind=8), allocatable:: hist(:,:) ! hist(bin, k) for each temperature k
  integer(kind=8):: hist_counts(1) ! dummy, reuse NT
  integer(kind=8):: nmeas_eq, nmeas_total ! number of measurements after thermalization and total
  integer:: ibin, Qi, i, rA, rB

  inputfile = '../inputs/part1.in'
  ! Suggested by co-pilot to allow input file as a command line argument
  call get_command_argument(1, inputfile)

  ! Extracts a short name from the input filename
  call make_tag_from_input(trim(inputfile), tag)

  ! Read input and initialize globals
  call read_input(trim(inputfile))
  call setr1279(seed)
  call init_globals()

  ! Initialize the lattice and bonds
  call alloc_lattice()
  call build_lattice()
  call alloc_bonds()
  call init_bonds()

  ! Initialize histogram
  allocate(hist(2*N+1, NT))
  hist(:,:) = 0_8
  nmeas_total = 0_8

  ! Allocate all the arrays for the spins, energies, maps, and statistics
  allocate(spinsA(N,NT), spinsB(N,NT))
  allocate(EA(NT), EB(NT))
  allocate(repA(NT), repB(NT), tempA(NT), tempB(NT))
  allocate(accA(NT-1), attA(NT-1), accB(NT-1), attB(NT-1))
  ! And post-transient swap statistics
  allocate(accA_eq(NT-1), attA_eq(NT-1), accB_eq(NT-1), attB_eq(NT-1))
  allocate(accA_old(NT-1), attA_old(NT-1), accB_old(NT-1), attB_old(NT-1))

  accA(:) = 0_8
  attA(:) = 0_8
  accB(:) = 0_8
  attB(:) = 0_8
  accA_eq(:) = 0_8
  attA_eq(:) = 0_8
  accB_eq(:) = 0_8
  attB_eq(:) = 0_8

  call init_family(spinsA, EA)
  call init_family(spinsB, EB)
  call init_maps(repA, tempA)
  call init_maps(repB, tempB)

  ! Output file names
  tsfile = 'timeseries_' // trim(tag) // '.dat'
  swfile = 'swap_stats_' // trim(tag) // '.dat'

  open(10, file=trim(tsfile), status='replace')
  write(10,'(A)') '# mcs  k  T(k)  Q(k)  Eavg(k)'

  ! This means the first exchange step will try pairs: (1,2), (3,4), (5,6), ...
  parity = 1

  ! Measure the swap acceptance after 10% of the total MCS: 
  mcs_therm = max(1000_8, nMCS / 10_8)
  print *, 'Post-transient swap statistics will be measured after MCS =', mcs_therm

  do mcs = 1_8, nMCS ! One pass of this loop corresponds to one Monte Carlo step for the whole PT simulation

    ! Metropolis sweeps for each replica at its current temperature
    do r = 1, NT
      call metropolis_sweep(spinsA(:,r), EA(r), beta_list(tempA(r)))
      call metropolis_sweep(spinsB(:,r), EB(r), beta_list(tempB(r)))
    enddo

    ! Every nsw steps, attempt exchanges between neighboring replicas in the family
    if (mod(mcs, nsw) == 0_8) then

      ! Save counters before the exchange step
      accA_old = accA
      attA_old = attA
      accB_old = accB
      attB_old = attB

      call attempt_exchange_family(EA, repA, tempA, parity, accA, attA)
      call attempt_exchange_family(EB, repB, tempB, parity, accB, attB)

      ! Accumulate only post-transient swap statistics
      if (mcs > mcs_therm) then
        accA_eq = accA_eq + (accA - accA_old)
        attA_eq = attA_eq + (attA - attA_old)
        accB_eq = accB_eq + (accB - accB_old)
        attB_eq = attB_eq + (attB - attB_old)
      endif

      ! The parity is flipped, so the next swap step uses the opposite set of neighboring pairs
      if (parity == 1) then
        parity = 2
      else
        parity = 1
      endif
    endif

    ! NOW: measure the overlap histogram for each temp
    if ((mod(mcs, nmeas) == 0_8).and.(mcs > mcs_therm)) then
      do k = 1, NT
        rA = repA(k)
        rB = repB(k)
        Qi = 0
        do i = 1, N
          Qi = Qi + spinsA(i,rA) * spinsB(i,rB)
        enddo
        ! Qi from -N to +N. so map to bin index 1..2N+1
        ibin = Qi + N + 1
        hist(ibin, k) = hist(ibin, k) + 1_8
      enddo
      nmeas_total = nmeas_total + 1_8
    endif

  enddo

  close(10)
  call write_swap_stats(swfile, accA, attA, accB, attB)

  ! Print post-transient acceptance rates on screen
  !print *
  print *, 'Swap acceptance rates after transient:'
  print *, ' pair      Tk          Tk+1        rateA       rateB'
  ok_all = .true.

  do k = 1, NT - 1
    rateA = safe_rate(accA_eq(k), attA_eq(k))
    rateB = safe_rate(accB_eq(k), attB_eq(k))

    write(*,'(I4,2X,F10.5,2X,F10.5,2X,F8.4,2X,F8.4)') &
         k, temp_list(k), temp_list(k+1), rateA, rateB

    if (rateA < 0.3d0 .or. rateB < 0.3d0) ok_all = .false.
  enddo

  print *
  if (ok_all) then
    print *, 'OK: all post-transient swap acceptance rates are >= 0.3'
  else
    print *, 'WARNING: some post-transient swap acceptance rates are < 0.3'
    print *, 'Try increasing NT (or nMCS)'
  endif

  deallocate(spinsA, spinsB, EA, EB, repA, repB, tempA, tempB, accA, attA, accB, attB)
  deallocate(accA_eq, attA_eq, accB_eq, attB_eq)
  deallocate(accA_old, attA_old, accB_old, attB_old)

  call dealloc_bonds()
  call dealloc_lattice()

contains

! Extract a short tag from the input filename for labeling output files
! suggested by co-pilot
  subroutine make_tag_from_input(fname, tag)
    implicit none
    character(*), intent(in)  :: fname
    character(*), intent(out) :: tag
    integer :: i1, i2

    i1 = index(fname, '/', back=.true.)
    if (i1 == 0) i1 = 0
    i2 = index(fname, '.', back=.true.)
    if (i2 == 0 .or. i2 <= i1) then
      tag = adjustl(fname(i1+1:))
    else
      tag = adjustl(fname(i1+1:i2-1))
    end if
  end subroutine make_tag_from_input

end program main