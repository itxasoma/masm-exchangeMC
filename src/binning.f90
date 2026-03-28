! Re-used from https://github.com/itxasoma/momo-McMC-2DIsing/blob/main/binning.f90 
module binning
  implicit none
contains
  subroutine average_block(X, block_size, av_X, sigma_X)
    implicit none
    double precision, intent(in):: X(:)
    integer, intent(in):: block_size
    double precision, intent(out):: av_X, sigma_X
    integer:: n_data, n_blocks
    integer:: i, beg_idx, end_idx
    double precision:: co
    double precision, allocatable:: block_X(:)

    ! Initialize variables
    n_data = size(X)
    n_blocks = n_data / block_size

    ! Integer division of how many blocks we have
    if (n_blocks <= 1) then
      av_X = sum(X) / dble(n_data)
      sigma_X = 0.d0
      return
    endif

    allocate(block_X(n_blocks))

    block_X(:) = 0.d0
    co = 0.d0
    av_X = 0.d0
    sigma_X = 0.d0

    !  Compute each block average (for every data-block)
    do i = 1, n_blocks
      beg_idx = 1 + (i - 1) * block_size
      end_idx = i * block_size
      block_X(i) = sum(X(beg_idx:end_idx)) / dble(block_size)
    enddo

    ! Average and error calculation
    av_X = sum(block_X) / dble(n_blocks)

    do i = 1, n_blocks
      co = co + (block_X(i) - av_X)**2.d0
    enddo

    co = co / (dble(n_blocks) * dble(n_blocks - 1))
    sigma_X = sqrt(co)

    deallocate(block_X)
  end subroutine average_block

end module binning