! CSR sparse matrix-vector product y := A*x written with newer language
! facilities.  Each variant changes ONE thing relative to the F77 baseline
! (spmv_f77.f) and keeps the arithmetic identical: per row, a sum starting
! from 0 of val(k)*x(col(k)) in storage order.  A perfect compiler would
! emit the same code for all of them.  The tag in brackets is the facility
! under test; report.py groups results by it.
module spmv_modern
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)

  !> CSR matrix as a derived type with allocatable components.
  type, public :: csr_matrix
    integer :: n = 0
    integer, allocatable :: rowptr(:), col(:)
    real(dp), allocatable :: val(:)
  contains
    procedure :: matvec => csrmv_tbp
    procedure, non_overridable :: row_dot => csr_row_dot
  end type

  public :: csrmv_ctrl, csrmv_ashape, csrmv_dotprod, csrmv_sum, csrmv_assoc
  public :: csrmv_dtype, csrmv_rowfn, csrmv_dc

contains

  !> [control] the F77 baseline transcribed unchanged (explicit-shape, DO
  !> loops).  Should run at 1.00; its spread is the noise floor.
  subroutine csrmv_ctrl(n, rowptr, col, val, x, y)
    integer, intent(in) :: n, rowptr(n+1), col(*)
    real(dp), intent(in) :: val(*), x(*)
    real(dp), intent(out) :: y(n)
    integer :: i, k
    real(dp) :: t
    do i = 1, n
      t = 0
      do k = rowptr(i), rowptr(i+1) - 1
        t = t + val(k)*x(col(k))
      end do
      y(i) = t
    end do
  end subroutine

  !> [assumed shape] baseline loops, all arrays assumed-shape dummies
  subroutine csrmv_ashape(rowptr, col, val, x, y)
    integer, intent(in) :: rowptr(:), col(:)
    real(dp), intent(in) :: val(:), x(:)
    real(dp), intent(out) :: y(:)
    integer :: i, k
    real(dp) :: t
    do i = 1, size(y)
      t = 0
      do k = rowptr(i), rowptr(i+1) - 1
        t = t + val(k)*x(col(k))
      end do
      y(i) = t
    end do
  end subroutine

  !> [intrinsic + vector subscript] row as dot_product with a gathered x
  subroutine csrmv_dotprod(n, rowptr, col, val, x, y)
    integer, intent(in) :: n, rowptr(n+1), col(*)
    real(dp), intent(in) :: val(*), x(*)
    real(dp), intent(out) :: y(n)
    integer :: i
    do i = 1, n
      y(i) = dot_product(val(rowptr(i):rowptr(i+1)-1), x(col(rowptr(i):rowptr(i+1)-1)))
    end do
  end subroutine

  !> [intrinsic + vector subscript] row as sum of an elemental product
  subroutine csrmv_sum(n, rowptr, col, val, x, y)
    integer, intent(in) :: n, rowptr(n+1), col(*)
    real(dp), intent(in) :: val(*), x(*)
    real(dp), intent(out) :: y(n)
    integer :: i
    do i = 1, n
      y(i) = sum(val(rowptr(i):rowptr(i+1)-1)*x(col(rowptr(i):rowptr(i+1)-1)))
    end do
  end subroutine

  !> [associate] row slices named with associate, loop over the slices
  subroutine csrmv_assoc(n, rowptr, col, val, x, y)
    integer, intent(in) :: n, rowptr(n+1), col(*)
    real(dp), intent(in) :: val(*), x(*)
    real(dp), intent(out) :: y(n)
    integer :: i, k
    real(dp) :: t
    do i = 1, n
      associate (c => col(rowptr(i):rowptr(i+1)-1), v => val(rowptr(i):rowptr(i+1)-1))
        t = 0
        do k = 1, size(v)
          t = t + v(k)*x(c(k))
        end do
      end associate
      y(i) = t
    end do
  end subroutine

  !> [derived type] baseline loops over allocatable components A%...
  subroutine csrmv_dtype(a, x, y)
    type(csr_matrix), intent(in) :: a
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)
    integer :: i, k
    real(dp) :: t
    do i = 1, a%n
      t = 0
      do k = a%rowptr(i), a%rowptr(i+1) - 1
        t = t + a%val(k)*x(a%col(k))
      end do
      y(i) = t
    end do
  end subroutine

  !> [type-bound] same body with a class(csr_matrix) passed object,
  !> called as A%matvec(x, y) on a non-polymorphic type(csr_matrix) object
  subroutine csrmv_tbp(this, x, y)
    class(csr_matrix), intent(in) :: this
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)
    integer :: i, k
    real(dp) :: t
    do i = 1, this%n
      t = 0
      do k = this%rowptr(i), this%rowptr(i+1) - 1
        t = t + this%val(k)*x(this%col(k))
      end do
      y(i) = t
    end do
  end subroutine

  !> [pure function] per-row dot product through a pure type-bound
  !> function in the same module (needs inlining)
  subroutine csrmv_rowfn(a, x, y)
    type(csr_matrix), intent(in) :: a
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)
    integer :: i
    do i = 1, a%n
      y(i) = a%row_dot(i, x)
    end do
  end subroutine

  pure real(dp) function csr_row_dot(this, i, x) result(t)
    class(csr_matrix), intent(in) :: this
    integer, intent(in) :: i
    real(dp), intent(in) :: x(*)
    integer :: k
    t = 0
    do k = this%rowptr(i), this%rowptr(i+1) - 1
      t = t + this%val(k)*x(this%col(k))
    end do
  end function

  !> [do concurrent] rows as a do concurrent with a block-local accumulator
  subroutine csrmv_dc(n, rowptr, col, val, x, y)
    integer, intent(in) :: n, rowptr(n+1), col(*)
    real(dp), intent(in) :: val(*), x(*)
    real(dp), intent(out) :: y(n)
    integer :: i
    do concurrent (i = 1:n)
      block
        integer :: k
        real(dp) :: t
        t = 0
        do k = rowptr(i), rowptr(i+1) - 1
          t = t + val(k)*x(col(k))
        end do
        y(i) = t
      end block
    end do
  end subroutine

end module spmv_modern
