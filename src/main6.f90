program main6
  use rng_wrapper
  use GLOBAL
  use lattice
  use bonds
  use exchangeMC
  implicit none

  integer, allocatable:: spinsA(:), spinsB(:)
  integer(kind=8):: mcs
  double precision:: EA, EB
  character(len=200):: inputfile, tsfile, tag

  inputfile = '../inputs/part6.in'
  call get_command_argument(1, inputfile)

  call make_tag_from_input(trim(inputfile), tag)

  call read_input(trim(inputfile))
  call setr1279(seed)
  call init_globals()
  call alloc_lattice()
  call build_lattice()
  call alloc_bonds()
  call init_bonds()

  allocate(spinsA(N), spinsB(N))
  call random_spins(spinsA)
  call total_energy(spinsA, EA)
  call random_spins(spinsB)
  call total_energy(spinsB, EB)

  tsfile = 'timeseries_' // trim(tag) // '.dat'

  open(10, file=trim(tsfile), status='replace')
  write(10,'(A)') '# mcs  k  T(k)  Q(k)  Eavg(k)'

  ! Main MC loop, for NT = 1 temperature
  do mcs = 1_8, nMCS
    call metropolis_sweep(spinsA, EA, beta_list(1))
    call metropolis_sweep(spinsB, EB, beta_list(1))
    if (mod(mcs, nmeas) == 0_8) then
      call measure_one_temp(mcs, 10, spinsA, spinsB, EA, EB)
    endif
  enddo

  close(10)

  deallocate(spinsA, spinsB)

  call dealloc_bonds()
  call dealloc_lattice()

contains

  subroutine measure_one_temp(mcs, unit_ts, spinsA, spinsB, EA, EB)
    implicit none
    integer(kind=8), intent(in):: mcs
    integer, intent(in):: unit_ts
    integer, intent(in):: spinsA(:), spinsB(:)
    double precision, intent(in):: EA, EB
    integer:: i, Q
    double precision:: Emean
    Q = 0
    do i = 1, N
      Q = Q + spinsA(i) * spinsB(i)
    enddo
    Emean = 0.5d0 * (EA + EB)
    write(unit_ts,'(I12,1X,I4,1X,F12.6,1X,I12,1X,F20.10)') mcs, 1, temp_list(1), Q, Emean
  end subroutine measure_one_temp

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

end program main6