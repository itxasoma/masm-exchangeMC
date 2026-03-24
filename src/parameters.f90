module GLOBAL 
  implicit none 
  integer:: d, L, N, z, NT, nsw, seed 
  integer(kind=8):: nMCS, nmeas 
  double precision:: Tmin, Tmax 
  double precision, allocatable:: temp_list(:), beta_list(:) 
 character(len=64):: bond_mode ! Nearest or all
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
    close(10) 
   bond_mode = adjustl(trim(bond_mode)) 
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
    allocate(temp_list(NT), beta_list(NT)) 
    if (NT == 1) then 
      temp_list(1) = Tmin 
    else 
      do k = 1, NT 
        temp_list(k) = Tmin + (Tmax - Tmin)*dble(k-1)/dble(NT-1) 
      enddo 
    endif 
    do k = 1, NT 
      beta_list(k) = 1.d0/temp_list(k) 
    enddo 

    print*, 'temp_list=', temp_list
    print*, 'beta_list=', beta_list

  end subroutine init_globals 

end module GLOBAL 
