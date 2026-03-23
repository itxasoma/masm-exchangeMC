module GLOBAL 
  implicit none 
  integer:: d, L, N, z, NT, nsw, seed ! Simulation parameters
  integer(kind=8):: nMCS, nmeas ! MC steps and MC measurement step number
  double precision:: Tmin, Tmax ! Temperature range
  double precision, allocatable:: temp_list(:), beta_list(:) ! Lists of T and betas
! character(len=64):: bond_mode ! Nearest or all
contains

! Read input parameters from a file
  subroutine read_input(inputfile)
    implicit none
    character(*), intent(in):: inputfile
    open(10, file=trim(inputfile), status='old') 
    read(10,*) d 
    read(10,*) L 
    read(10,*) Tmin 
    read(10,*) Tmax 
    read(10,*) NT 
    read(10,*) nsw 
    read(10,*) nmeas 
    read(10,*) nMCS 
    read(10,*) seed 
!    read(10,'(A)') bond_mode ! New
    close(10) 
!    bond_mode = adjustl(trim(bond_mode)) ! New
  end subroutine read_input

! Initialize global variables and allocate T, beta lists
  subroutine init_globals() 
    implicit none 
    integer:: k, idir

    !N = 1 
    !do idir = 1, d 
    !  N = N*L 
    !enddo 

    !z = 2*d 

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
  end subroutine init_globals 

end module GLOBAL 
