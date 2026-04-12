module GLOBAL 
  implicit none 
  integer:: d, L, N, z, NT, nsw, seed, ios 
  integer(kind=8):: nMCS, nmeas 
  double precision:: Tmin, Tmax 
  double precision, allocatable:: temp_list(:), beta_list(:) 
  character(len=64):: bond_mode ! Nearest or all
  character(len=32):: temp_mode ! Linear spacing or other for temps
  character(len=200):: temp_file  ! File with list of temps
contains

! Read input parameters from a file
  subroutine read_input(inputfile)
    implicit none
    character(*), intent(in):: inputfile
    open(10, file=trim(inputfile), status='old') 
    read(10,*) d                    ! space dimension
    read(10,*) L                    ! linear size of the lattice     
    read(10,*) Tmin                 ! min temp
    read(10,*) Tmax                 ! max temp
    read(10,*) NT                   ! num of replicas (1 for each T)
    read(10,*) nsw                  ! num of MC steps between swap attempts
    read(10,*) nmeas                ! num of measurement steps
    read(10,*) nMCS                 ! num of total MCsteps 
    read(10,*) seed                 ! rng seed
    read(10,'(A)') bond_mode        ! bond type: gaussian or ferro
    bond_mode = adjustl(trim(bond_mode)) 

    ! VS copilot suggestion: read temp_mode, if not given, keep default 'linear'
    temp_mode = 'linear'
    temp_file = ''

    read(10,'(A)', iostat=ios) temp_mode 
    if (ios == 0) then
      temp_mode = adjustl(trim(temp_mode))
      if (trim(temp_mode) == 'file') then
        read(10,'(A)', iostat=ios) temp_file
        if (ios /= 0) stop 'temp_mode=file but no temperature file given'
        temp_file = adjustl(trim(temp_file))
      endif
    endif

    close(10) 
  end subroutine read_input

! Initialize global variables and allocate T, beta lists
  subroutine init_globals() 
    implicit none 
    integer:: k

    ! Total num of spins in the lattice:
    N = L**d
    print*, 'N=', N

    ! Coordination number:
    ! Num of nn for square/cubic lattices with pbcs
    z = 2*d 
    print*, 'z=', z

    ! The arrays' dimensions are NT: number of replicas.
    ! Tk are linearly spaced between Tmin and Tmax.
    if (allocated(temp_list)) deallocate(temp_list)
    if (allocated(beta_list)) deallocate(beta_list)
    allocate(temp_list(NT), beta_list(NT))

    if (trim(temp_mode) == 'file') then
      call read_temperature_file(trim(temp_file))
    else
      if (NT == 1) then
        temp_list(1) = Tmin
      else
        do k = 1, NT
          temp_list(k) = Tmin + (Tmax - Tmin) * dble(k-1) / dble(NT-1)
        enddo
      endif
    endif

    do k = 1, NT
      beta_list(k) = 1.d0 / temp_list(k)
    enddo

    print*, 'N=', N
    print*, 'z=', z
    print*, 'temp_list=', temp_list

  end subroutine init_globals 
  
  subroutine read_temperature_file(filename)
    implicit none
    character(*), intent(in):: filename
    integer:: k
    open(11, file=trim(filename), status='old')
    do k = 1, NT
      read(11,*) temp_list(k)
    enddo
    close(11)
  end subroutine read_temperature_file

end module GLOBAL 
