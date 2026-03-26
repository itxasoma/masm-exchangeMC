! Re-used from https://github.com/itxasoma/momo-McMC-2DIsing/blob/main/binning.f90 
module binning
  implicit none
contains
  subroutine average_block(X, block_size, av_X, sigma_X)
    implicit none
    double precision, intent(in):: X(:)
    integer, intent(in):: block_size
    double precision, intent(out):: av_X, sigma_X
    integer:: n_data, n_blocks, i, start_idx
    double precision, allocatable:: block_X(:)

    n_data = size(X)
    n_blocks = n_data / block_size

    if (n_blocks <= 1) then
      av_X = sum(X) / dble(n_data)
      sigma_X = 0.d0
      return
    end if

    allocate(block_X(n_blocks))

    do i = 1, n_blocks
      start_idx = 1 + (i-1) * block_size
      block_X(i) = sum(X(start_idx:start_idx + block_size - 1)) / dble(block_size)
    end do

    av_X = sum(block_X) / dble(n_blocks)
    sigma_X = sqrt(max(0.d0, sum((block_X - av_X)**2) / dble(n_blocks * (n_blocks - 1))))

    deallocate(block_X)
  end subroutine average_block

end module binning