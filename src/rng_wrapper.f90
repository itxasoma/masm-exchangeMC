! To call the external RNG functions cleanly
module rng_wrapper
  implicit none
  interface
    subroutine setr1279(iseed)
      integer, intent(in):: iseed
    end subroutine setr1279

    real function r1279()
    end function r1279

    integer function ir1279()
    end function ir1279

    integer function ir1279range(imin, imax)
      integer, intent(in):: imin, imax
    end function ir1279range
  end interface
end module rng_wrapper
