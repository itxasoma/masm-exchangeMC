program main
  use rng_wrapper
  use GLOBAL
  use lattice
  use bonds
  use exchangeMC
  implicit none

  integer:: r, parity, i, k, ibin, rA, rB, Qi
  integer, allocatable:: spinsA(:,:), spinsB(:,:)
  integer, allocatable:: repA(:), repB(:), tempA(:), tempB(:)
  integer(kind=8):: mcs, nmeas_total
  integer(kind=8), allocatable:: accA(:), attA(:), accB(:), attB(:)
  
  integer(kind=8), allocatable:: accA_eq(:), attA_eq(:), accB_eq(:), attB_eq(:)
  integer(kind=8), allocatable:: accA_old(:), attA_old(:), accB_old(:), attB_old(:)
  integer(kind=8), allocatable:: hist(:,:)   ! hist(ibin, k)
  double precision, allocatable:: EA(:), EB(:)
  character(len=200):: inputfile, histfile, swfile, tag
  integer(kind=8):: mcs_therm
  double precision:: rateA, rateB
  logical:: ok_all

  inputfile = '../inputs/part1.in'
  call get_command_argument(1, inputfile)
  call make_tag_from_input(trim(inputfile), tag)

  call read_input(trim(inputfile))
  call setr1279(seed)
  call init_globals()

  call alloc_lattice()
  call build_lattice()
  call alloc_bonds()
  call init_bonds()

  allocate(spinsA(N,NT), spinsB(N,NT))
  allocate(EA(NT), EB(NT))
  allocate(repA(NT), repB(NT), tempA(NT), tempB(NT))
  allocate(accA(NT-1), attA(NT-1), accB(NT-1), attB(NT-1))
  allocate(accA_eq(NT-1), attA_eq(NT-1), accB_eq(NT-1), attB_eq(NT-1))
  allocate(accA_old(NT-1), attA_old(NT-1), accB_old(NT-1), attB_old(NT-1))
  ! Histogram: bins for Q = -N .. +N, i.e. 2N+1 bins, for each temperature slot
  allocate(hist(2*N+1, NT))

  accA(:) = 0_8;  attA(:) = 0_8
  accB(:) = 0_8;  attB(:) = 0_8
  accA_eq(:) = 0_8;  attA_eq(:) = 0_8
  accB_eq(:) = 0_8;  attB_eq(:) = 0_8
  hist(:,:) = 0_8
  nmeas_total = 0_8

  call init_family(spinsA, EA)
  call init_family(spinsB, EB)
  call init_maps(repA, tempA)
  call init_maps(repB, tempB)

  histfile = 'histogram_' // trim(tag) // '.dat'
  swfile   = 'swap_stats_' // trim(tag) // '.dat'

  parity    = 1
  mcs_therm = max(3000_8, nMCS / 10_8)   ! discard first 3000 MCS (assignment requirement)
  print *, 'Thermalization cut at MCS =', mcs_therm

  do mcs = 1_8, nMCS

    do r = 1, NT
      call metropolis_sweep(spinsA(:,r), EA(r), beta_list(tempA(r)))
      call metropolis_sweep(spinsB(:,r), EB(r), beta_list(tempB(r)))
    enddo

    if (mod(mcs, nsw) == 0_8) then
      accA_old = accA;  attA_old = attA
      accB_old = accB;  attB_old = attB

      call attempt_exchange_family(EA, repA, tempA, parity, accA, attA)
      call attempt_exchange_family(EB, repB, tempB, parity, accB, attB)

      if (mcs > mcs_therm) then
        accA_eq = accA_eq + (accA - accA_old)
        attA_eq = attA_eq + (attA - attA_old)
        accB_eq = accB_eq + (accB - accB_old)
        attB_eq = attB_eq + (attB - attB_old)
      endif

      parity = 3 - parity   ! toggles between 1 and 2
    endif

    ! Accumulate histogram only after thermalization
    if (mod(mcs, nmeas) == 0_8 .and. mcs > mcs_therm) then
      do k = 1, NT
        rA = repA(k)
        rB = repB(k)
        Qi = 0
        do i = 1, N
          Qi = Qi + spinsA(i,rA) * spinsB(i,rB)
        end do
        ibin = Qi + N + 1          ! maps Q in [-N,N] -> index in [1, 2N+1]
        hist(ibin, k) = hist(ibin, k) + 1_8
      end do
      nmeas_total = nmeas_total + 1_8
    endif

  enddo

  ! Write compact histogram file (NT * (2N+1) lines ~ a few KB)
  open(10, file=trim(histfile), status='replace')
  write(10,'(A,I12)') '# nmeas_total = ', nmeas_total
  write(10,'(A)') '# k  T  Q  count'
  do k = 1, NT
    do ibin = 1, 2*N+1
      Qi = ibin - N - 1
      write(10,'(I4,1X,F12.6,1X,I6,1X,I12)') k, temp_list(k), Qi, hist(ibin,k)
    end do
  end do
  close(10)

  call write_swap_stats(swfile, accA, attA, accB, attB)

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
  endif

  deallocate(spinsA, spinsB, EA, EB, repA, repB, tempA, tempB)
  deallocate(accA, attA, accB, attB)
  deallocate(accA_eq, attA_eq, accB_eq, attB_eq)
  deallocate(accA_old, attA_old, accB_old, attB_old)
  deallocate(hist)
  call dealloc_bonds()
  call dealloc_lattice()

contains

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