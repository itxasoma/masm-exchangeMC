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

  ! Allocate all the arrays for the spins, energies, maps, and statistics
  allocate(spinsA(N,NT), spinsB(N,NT))
  allocate(EA(NT), EB(NT))
  allocate(repA(NT), repB(NT), tempA(NT), tempB(NT))
  allocate(accA(NT-1), attA(NT-1), accB(NT-1), attB(NT-1))

  accA(:) = 0_8
  attA(:) = 0_8
  accB(:) = 0_8
  attB(:) = 0_8

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

  do mcs = 1_8, nMCS ! One pass of this loop corresponds to one Monte Carlo step for the whole PT simulation

    ! Metropolis sweeps for each replica at its current temperature
    do r = 1, NT
      call metropolis_sweep(spinsA(:,r), EA(r), beta_list(tempA(r)))
      call metropolis_sweep(spinsB(:,r), EB(r), beta_list(tempB(r)))
    enddo

    ! Every nsw steps, attempt exchanges between neighboring replicas in the family
    if (mod(mcs, nsw) == 0_8) then
      call attempt_exchange_family(EA, repA, tempA, parity, accA, attA)
      call attempt_exchange_family(EB, repB, tempB, parity, accB, attB)
      ! The parity is flipped, so the next swap step uses the opposite set of neighboring pairs
      if (parity == 1) then
        parity = 2
      else
        parity = 1
      endif
    endif

    ! Write observables every nmeas steps
    if (mod(mcs, nmeas) == 0_8) then
      call measure_all(mcs, 10, spinsA, spinsB, EA, EB, repA, repB)
    endif

  enddo

  close(10)
  call write_swap_stats(swfile, accA, attA, accB, attB)

  deallocate(spinsA, spinsB, EA, EB, repA, repB, tempA, tempB, accA, attA, accB, attB)
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
