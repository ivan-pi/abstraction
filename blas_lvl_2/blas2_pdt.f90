! [PDT] GEMV N on a parameterized derived type whose length parameters
! are the storage bounds.  Kept in its own file because PDT support is
! uneven; run.sh builds without it if this file does not compile.
module blas2_pdt
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)

  type, public :: pmat(ld, nc)
    integer, len :: ld, nc
    integer :: m = 0
    real(dp) :: a(ld, nc)
  end type

  public :: gemvn_pdt

contains

  subroutine gemvn_pdt(mat, x, y)
    type(pmat(*,*)), intent(in) :: mat
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)
    integer :: i, j
    real(dp) :: temp
    do i = 1, mat%m
      y(i) = 0
    end do
    do j = 1, mat%nc
      temp = x(j)
      do i = 1, mat%m
        y(i) = y(i) + temp*mat%a(i, j)
      end do
    end do
  end subroutine

end module blas2_pdt
