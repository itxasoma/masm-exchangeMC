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
      ! Metropolis: 
      if (dE <= 0.d0) then
        spin(ispin) = -spin(ispin)
        E = E + dE
      else if (dble(r1279()) < exp(-beta_replica * dE)) then
        spin(ispin) = -spin(ispin)
        E = E + dE
      endif

    enddo
  end subroutine metropolis_sweep

! Initialize the family of replicas with random spins and compute their energies
  subroutine init_family(spins, Erep)
    implicit none
    integer, intent(out):: spins(N,NT)
    double precision, intent(out):: Erep(NT)
    integer:: r
    do r = 1, NT
      call random_spins(spins(:,r))
      call total_energy(spins(:,r), Erep(r))
    enddo
  end subroutine init_family

! Initialize with replica r at temperature r
! It sets up the initial correspondence between replicas and temperature slots. 
! At the very beginning, nothing has been exchanged yet, so replica 1 is at temperature 1,
! replica 2 is at temperature 2, and so on.
  subroutine init_maps(rep_at_temp, temp_of_rep)
    implicit none
    integer, intent(out):: rep_at_temp(NT), temp_of_rep(NT)
    integer:: k
    do k = 1, NT
      rep_at_temp(k) = k ! Which replica is currently sitting at temperature slot k?
      temp_of_rep(k) = k ! At which temperature slot is replica r currently sitting?
    enddo
  end subroutine init_maps


! Attempt exchanges between neighboring replicas (k, k+1) in the family.
! Update the maps and acceptance statsics for each pair of neighboring temperatures.
  subroutine attempt_exchange_family(Erep, rep_at_temp, temp_of_rep, parity, acc, att)
    implicit none
    double precision, intent(in):: Erep(NT)
    integer, intent(inout):: rep_at_temp(NT), temp_of_rep(NT)
    integer, intent(in):: parity
    integer(kind=8), intent(inout):: acc(NT-1), att(NT-1)
    integer:: k, r1, r2, t1, t2, tmprep
    double precision:: expo

    ! Assignment’s requirement to alternate even and odd neighboring pairs
    do k = parity, NT - 1, 2 ! Attempts between pairs (1,2), (3,4), or (2,3), (4,5), etc
      r1 = rep_at_temp(k)
      r2 = rep_at_temp(k+1)

      t1 = temp_of_rep(r1)
      t2 = temp_of_rep(r2)

      ! a = min(1, exp((beta(t1) - beta(t2)) * (Erep(r1) - Erep(r2))))
      expo = (beta_list(t1) - beta_list(t2)) * (Erep(r1) - Erep(r2))

      ! attempt counter
      att(k) = att(k) + 1_8

      if (expo >= 0.d0) then
        tmprep = rep_at_temp(k)
        rep_at_temp(k) = rep_at_temp(k+1)
        rep_at_temp(k+1) = tmprep
        temp_of_rep(r1) = t2
        temp_of_rep(r2) = t1
        acc(k) = acc(k) + 1_8

      else if (dble(r1279()) < exp(expo)) then
        tmprep = rep_at_temp(k)
        rep_at_temp(k) = rep_at_temp(k+1)
        rep_at_temp(k+1) = tmprep
        temp_of_rep(r1) = t2
        temp_of_rep(r2) = t1
        acc(k) = acc(k) + 1_8
      endif

    enddo
  end subroutine attempt_exchange_family

! Measure overlap Q(k) and average energy E(k) by temperature slot k
  subroutine measure_all(mcs, unit_ts, spinsA, spinsB, EA, EB, repA, repB)
    implicit none
    integer(kind=8), intent(in):: mcs
    integer, intent(in):: unit_ts
    integer, intent(in):: spinsA(N,NT), spinsB(N,NT)
    double precision, intent(in):: EA(NT), EB(NT)
    integer, intent(in):: repA(NT), repB(NT)
    integer:: k, rA, rB, i, Q
    double precision:: Emean

    do k = 1, NT
      rA = repA(k)
      rB = repB(k)
      Q = 0
      do i = 1, N
        Q = Q + spinsA(i,rA) * spinsB(i,rB)
      end do
      Emean = 0.5d0 * (EA(rA) + EB(rB))
      write(unit_ts,'(I12,1X,I4,1X,F12.6,1X,I12,1X,F20.10)') mcs, k, temp_list(k), Q, Emean
    end do
  end subroutine measure_all

! A safe way to compute acceptance rates without risking division by zero
! Suggested by co-pilot
  double precision function safe_rate(acc, att)
    implicit none
    integer(kind=8), intent(in):: acc, att
    if (att <= 0_8) then
      safe_rate = 0.d0
    else
      safe_rate = dble(acc) / dble(att)
    end if
  end function safe_rate

! Write all the statistics about swap attempts
  subroutine write_swap_stats(filename, accA, attA, accB, attB)
    implicit none
    character(*), intent(in) :: filename
    integer(kind=8), intent(in) :: accA(NT-1), attA(NT-1), accB(NT-1), attB(NT-1)
    integer :: k, unit_sw

    unit_sw = 92
    open(unit_sw, file=trim(filename), status='replace')
    write(unit_sw,'(A)') '# pair Tk Tk+1 accA attA rateA accB attB rateB'
    do k = 1, NT - 1
      write(unit_sw,'(I4,1X,F10.5,1X,F10.5,1X,I12,1X,I12,1X,F10.5,1X,I12,1X,I12,1X,F10.5)') &
           k, temp_list(k), temp_list(k+1), accA(k), attA(k), safe_rate(accA(k),attA(k)), &
           accB(k), attB(k), safe_rate(accB(k),attB(k))
    end do
    close(unit_sw)
  end subroutine write_swap_stats


end module exchangeMC

