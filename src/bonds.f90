module bonds
  use GLOBAL
  use lattice
  use rng_wrapper
  implicit none
  double precision, allocatable:: Jij(:,:)
contains

! Allocate/deallocate the bond array
  subroutine alloc_bonds()
    implicit none
    allocate(Jij(N, z))
  end subroutine alloc_bonds

  subroutine dealloc_bonds()
    implicit none
    if (allocated(Jij)) deallocate(Jij)
  end subroutine dealloc_bonds

! Initialize bond array with gaussian or ferro
  subroutine init_bonds()
    implicit none
    Jij(:,:) = 0.d0
    select case (trim(adjustl(bond_mode)))
    case ('gaussian')
      call set_gaussian_bonds()
    case ('ferro')
      call set_ferro_bonds()
    case default
      stop 'bond_mode must be gaussian or ferro'
    end select
  end subroutine init_bonds

  subroutine set_ferro_bonds()
    implicit none
    integer:: i
    do i = 1, N
      ! Only positive directions: +x=1, +y=3, +z=5
      call set_bond_pair(i, 1, 2, 1.d0)
      if (d >= 2) call set_bond_pair(i, 3, 4, 1.d0)
      if (d >= 3) call set_bond_pair(i, 5, 6, 1.d0)
    end do
  end subroutine set_ferro_bonds

  subroutine set_gaussian_bonds()
    implicit none
    integer :: i
    double precision :: val
    do i = 1, N
      ! +x direction: only set if nbr(i,1) > i to avoid double-setting
      if (nbr(i,1) > i .or. nbr(i,1) < i) then
        ! simpler: only set when i < j (forward bond)
      endif
    end do
  end subroutine set_gaussian_bonds

! Set the bond value for a pair of neighbors i and j in the given direction
  subroutine set_bond_pair(i, dir_pos, dir_neg, val)
    implicit none
    integer, intent(in):: i, dir_pos, dir_neg
    double precision, intent(in):: val
    integer:: j
    j = nbr(i, dir_pos)
    Jij(i, dir_pos) = val
    Jij(j, dir_neg) = val
  end subroutine set_bond_pair

  double precision function gauss_rand()
    implicit none
    double precision:: u1, u2, v1, v2, s
10  continue
    u1 = dble(r1279())
    u2 = dble(r1279())
    v1 = 2.d0*u1 - 1.d0
    v2 = 2.d0*u2 - 1.d0
    s = v1*v1 + v2*v2
    if (s >= 1.d0 .or. s <= 1.d-14) goto 10
    gauss_rand = v1 * sqrt(-2.d0*log(s)/s)
  end function gauss_rand

end module bonds
