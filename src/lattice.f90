module lattice
  use GLOBAL
  implicit none
  integer, allocatable:: nbr(:,:)   ! nbr(i,k): k-th neighbour of site i
  integer, allocatable:: ina(:,:)   ! periodic shift array
contains

! Allocate/deallocate the neighbor array
  subroutine alloc_lattice()
    implicit none
    allocate(nbr(N, z))
    allocate(ina(0:1, L))
  end subroutine

  subroutine dealloc_lattice()
    implicit none
    if (allocated(nbr)) deallocate(nbr)
    if (allocated(ina)) deallocate(ina)
  end subroutine

! Apply pbcs
  subroutine pbc()
    implicit none
    integer:: i
      do i = 1, L
          ina(0,i) = i-1
          ina(1,i) = i+1
      enddo
      ina(0,1) = L
      ina(1,L) = 1
  end subroutine

! Build neighbor lists for 2d and 3d lattices + pbc 
  subroutine build_lattice()
    implicit none
    call pbc()
    if (d == 2) then
      call square_lattice()
    elseif (d == 3) then
      call cubic_lattice()
    else
      stop 'Choose d=2 or d=3.'
    endif
  end subroutine

! 2d square lattice + pbc
  subroutine square_lattice()
    implicit none
    integer:: i
    integer:: x, y
    i = 0
      do y = 1, L
        do x = 1, L
          i = i + 1
            nbr(i,1) = ina(1,x) + L*(y-1)
            nbr(i,2) = ina(0,x) + L*(y-1)
            nbr(i,3) = x + L*(ina(1,y)-1)
            nbr(i,4) = x + L*(ina(0,y)-1)
        enddo
      enddo
  end subroutine

! 3d cubic lattice + pbc
  subroutine cubic_lattice()
      implicit none
      integer:: i
      integer:: x, y, iz
      i = 0
      do iz = 1, L
          do y = 1, L
              do x = 1, L
                  i = i + 1
                  nbr(i,1) = idx3(ina(1,x), y, iz)
                  nbr(i,2) = idx3(ina(0,x), y, iz)
                  nbr(i,3) = idx3(x, ina(1,y), iz)
                  nbr(i,4) = idx3(x, ina(0,y), iz)
                  nbr(i,5) = idx3(x, y, ina(1,iz))
                  nbr(i,6) = idx3(x, y, ina(0,iz))
              enddo
          enddo
      enddo
  end subroutine

  integer function idx3(x, y, iz)
    implicit none
    integer, intent(in):: x, y, iz
    idx3 = x + L*(y-1) + L*L*(iz-1)
  end function

end module
