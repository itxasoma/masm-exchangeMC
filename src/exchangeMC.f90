module exchangeMC
  use GLOBAL
  use lattice
  use bonds
  use rng_wrapper
  implicit none
contains

! Initialize the spins randomly and compute initial E for each replica
  subroutine random_spins(spin)
    implicit none
    integer, intent(out):: spin(:)
    integer:: i
    do i = 1, size(spin)
      if (r1279() < 0.5) then
        spin(i) = 1
      else
        spin(i) = -1
      endif
    enddo
  end subroutine random_spins

  subroutine total_energy(spin, E)
    implicit none
    integer, intent(in):: spin(:)
    double precision, intent(out):: E
    integer:: i, k
    E = 0.d0
    do i = 1, N
      do k = 1, z
        E = E - 0.5d0 * Jij(i,k) * dble(spin(i) * spin(nbr(i,k)))
      enddo
    enddo
  end subroutine total_energy

! Energy change if we flip spin ispin 
  double precision function delta_energy(spin, ispin)
    implicit none
    integer, intent(in):: spin(:), ispin
    integer:: k
    double precision:: h
    h = 0.d0
    do k = 1, z
      h = h + Jij(ispin,k) * dble(spin(nbr(ispin,k)))
    end do
    delta_energy = 2.d0 * dble(spin(ispin)) * h
  end function delta_energy

! Perform exchange MC sweeps for all replicas
! N attempts per sweep, where N is the total num of spins in the lattice
  subroutine metropolis_sweep(spin, E, beta_replica)
    implicit none
    integer, intent(inout):: spin(:)
    double precision, intent(inout):: E
    double precision, intent(in):: beta_replica
    integer:: attempt, ispin
    double precision:: dE

    do attempt = 1, N
      ispin = min(N, int(dble(N) * dble(r1279())) + 1)
      dE = delta_energy(spin, ispin)

      if (dE <= 0.d0) then
        spin(ispin) = -spin(ispin)
        E = E + dE
      else if (dble(r1279()) < exp(-beta_replica * dE)) then
        spin(ispin) = -spin(ispin)
        E = E + dE
      endif
    enddo
  end subroutine metropolis_sweep


end module exchangeMC

