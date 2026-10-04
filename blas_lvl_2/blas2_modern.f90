! BLAS level-2 kernels written with newer language facilities.
!
! Each variant changes ONE thing relative to the F77 baseline in
! blas2_f77.f and keeps the arithmetic (operand order, parenthesisation,
! summation order) identical, so a perfect compiler would emit the same
! code.  The tag in brackets is the facility under test; report.py groups
! results by it.
module blas2_modern
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)

  !> Matrix wrapper: shape plus an allocatable (lda x n) payload.
  type, public :: dmat
    integer :: m = 0, n = 0
    real(dp), allocatable :: a(:,:)
  contains
    procedure :: mv => gemvn_tbp
  end type

  !> Upper triangle in column-packed storage, element access via a
  !> type-bound accessor instead of hand-maintained packed indices.
  type, public :: packed_upper
    integer :: n = 0
    real(dp), allocatable :: ap(:)
  contains
    procedure, non_overridable :: at => packed_at
  end type

  public :: gemvn_arrsec, gemvn_ashape, gemvn_matmul, gemvn_dtype
  public :: gemvt_sum, gemvt_dotprod, gemvt_matmul_vm, gemvt_matmul_tr
  public :: ger_arrsec, ger_dc_ji, ger_dc_ij, ger_spread
  public :: trsv_arrsec, trsv_arrsec_tmp
  public :: symv_arrsec_sum, symv_dotprod
  public :: tpmv_arrsec, tpmv_accessor

contains

  ! ------------------------------------------------------------------ GEMV N

  !> [array sections] column axpy as a section update
  subroutine gemvn_arrsec(m, n, a, lda, x, y)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    integer :: j
    y(1:m) = 0
    do j = 1, n
      y(1:m) = y(1:m) + x(j)*a(1:m,j)
    end do
  end subroutine

  !> [assumed shape] same loops, but A, x, y are assumed-shape dummies and
  !> the caller passes the non-contiguous section a(1:m,1:n) of an lda array.
  subroutine gemvn_ashape(a, x, y)
    real(dp), intent(in) :: a(:,:), x(:)
    real(dp), intent(out) :: y(:)
    integer :: i, j
    real(dp) :: temp
    do i = 1, size(a, 1)
      y(i) = 0
    end do
    do j = 1, size(a, 2)
      temp = x(j)
      do i = 1, size(a, 1)
        y(i) = y(i) + temp*a(i, j)
      end do
    end do
  end subroutine

  !> [intrinsic] matmul on explicit-shape sections
  subroutine gemvn_matmul(m, n, a, lda, x, y)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    y(1:m) = matmul(a(1:m,1:n), x(1:n))
  end subroutine

  !> [derived type] loops over the allocatable component mat%a
  subroutine gemvn_dtype(mat, x, y)
    type(dmat), intent(in) :: mat
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)
    integer :: i, j
    real(dp) :: temp
    do i = 1, mat%m
      y(i) = 0
    end do
    do j = 1, mat%n
      temp = x(j)
      do i = 1, mat%m
        y(i) = y(i) + temp*mat%a(i, j)
      end do
    end do
  end subroutine

  !> [type-bound] same body with a polymorphic class(dmat) passed object,
  !> called as mat%mv(x, y) on a non-polymorphic type(dmat) object.
  subroutine gemvn_tbp(this, x, y)
    class(dmat), intent(in) :: this
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)
    integer :: i, j
    real(dp) :: temp
    do i = 1, this%m
      y(i) = 0
    end do
    do j = 1, this%n
      temp = x(j)
      do i = 1, this%m
        y(i) = y(i) + temp*this%a(i, j)
      end do
    end do
  end subroutine

  ! ------------------------------------------------------------------ GEMV T

  !> [intrinsic] column dot product as sum of an elemental product
  subroutine gemvt_sum(m, n, a, lda, x, y)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    integer :: j
    do j = 1, n
      y(j) = sum(a(1:m,j)*x(1:m))
    end do
  end subroutine

  !> [intrinsic] column dot product via dot_product
  subroutine gemvt_dotprod(m, n, a, lda, x, y)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    integer :: j
    do j = 1, n
      y(j) = dot_product(a(1:m,j), x(1:m))
    end do
  end subroutine

  !> [intrinsic] vector-matrix matmul: y = x^T A
  subroutine gemvt_matmul_vm(m, n, a, lda, x, y)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    y(1:n) = matmul(x(1:m), a(1:m,1:n))
  end subroutine

  !> [intrinsic] matmul(transpose(A), x): needs transpose folding
  subroutine gemvt_matmul_tr(m, n, a, lda, x, y)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    y(1:n) = matmul(transpose(a(1:m,1:n)), x(1:m))
  end subroutine

  ! --------------------------------------------------------------------- GER

  !> [array sections] column update as a section
  subroutine ger_arrsec(m, n, alpha, x, y, a, lda)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda,*)
    integer :: j
    real(dp) :: temp
    do j = 1, n
      temp = alpha*y(j)
      a(1:m,j) = a(1:m,j) + x(1:m)*temp
    end do
  end subroutine

  !> [do concurrent] index order (j, i): column-major friendly as written
  subroutine ger_dc_ji(m, n, alpha, x, y, a, lda)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda,*)
    integer :: i, j
    do concurrent (j = 1:n, i = 1:m)
      a(i,j) = a(i,j) + x(i)*(alpha*y(j))
    end do
  end subroutine

  !> [do concurrent] index order (i, j): iteration order is unspecified, so
  !> a good compiler should produce the same code as (j, i).
  subroutine ger_dc_ij(m, n, alpha, x, y, a, lda)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda,*)
    integer :: i, j
    do concurrent (i = 1:m, j = 1:n)
      a(i,j) = a(i,j) + x(i)*(alpha*y(j))
    end do
  end subroutine

  !> [intrinsic] outer product via two spread temporaries (removable)
  subroutine ger_spread(m, n, alpha, x, y, a, lda)
    integer, intent(in) :: m, n, lda
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda,*)
    a(1:m,1:n) = a(1:m,1:n) + spread(x(1:m), 2, n)*spread(alpha*y(1:n), 1, m)
  end subroutine

  ! -------------------------------------------------------------------- TRSV

  !> [array sections] x(j) appears on the RHS of an assignment to
  !> x(1:j-1); the compiler must prove no overlap to avoid a temporary.
  subroutine trsv_arrsec(n, a, lda, x)
    integer, intent(in) :: n, lda
    real(dp), intent(in) :: a(lda,*)
    real(dp), intent(inout) :: x(*)
    integer :: j
    do j = n, 1, -1
      x(j) = x(j)/a(j,j)
      x(1:j-1) = x(1:j-1) - x(j)*a(1:j-1,j)
    end do
  end subroutine

  !> [array sections] same, with x(j) hoisted into a scalar by hand
  subroutine trsv_arrsec_tmp(n, a, lda, x)
    integer, intent(in) :: n, lda
    real(dp), intent(in) :: a(lda,*)
    real(dp), intent(inout) :: x(*)
    integer :: j
    real(dp) :: temp
    do j = n, 1, -1
      x(j) = x(j)/a(j,j)
      temp = x(j)
      x(1:j-1) = x(1:j-1) - temp*a(1:j-1,j)
    end do
  end subroutine

  ! -------------------------------------------------------------------- SYMV

  !> [array sections] the fused axpy+dot loop split into two statements;
  !> recovering one pass over a(1:j-1,j) needs loop fusion.
  subroutine symv_arrsec_sum(n, a, lda, x, y)
    integer, intent(in) :: n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    integer :: j
    real(dp) :: temp1, temp2
    y(1:n) = 0
    do j = 1, n
      temp1 = x(j)
      y(1:j-1) = y(1:j-1) + temp1*a(1:j-1,j)
      temp2 = sum(a(1:j-1,j)*x(1:j-1))
      y(j) = y(j) + temp1*a(j,j) + temp2
    end do
  end subroutine

  !> [intrinsic] as above with dot_product
  subroutine symv_dotprod(n, a, lda, x, y)
    integer, intent(in) :: n, lda
    real(dp), intent(in) :: a(lda,*), x(*)
    real(dp), intent(out) :: y(*)
    integer :: j
    real(dp) :: temp1, temp2
    y(1:n) = 0
    do j = 1, n
      temp1 = x(j)
      y(1:j-1) = y(1:j-1) + temp1*a(1:j-1,j)
      temp2 = dot_product(a(1:j-1,j), x(1:j-1))
      y(j) = y(j) + temp1*a(j,j) + temp2
    end do
  end subroutine

  ! -------------------------------------------------------------------- TPMV

  !> [array sections] packed column as a section with computed bounds
  subroutine tpmv_arrsec(n, ap, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: ap(*), x(*)
    real(dp), intent(out) :: y(*)
    integer :: j, kk
    y(1:n) = 0
    kk = 1
    do j = 1, n
      y(1:j) = y(1:j) + x(j)*ap(kk:kk+j-1)
      kk = kk + j
    end do
  end subroutine

  !> [accessor] element access through a pure type-bound function that
  !> recomputes the packed index; needs inlining + strength reduction.
  subroutine tpmv_accessor(u, x, y)
    type(packed_upper), intent(in) :: u
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: y(*)
    integer :: i, j
    real(dp) :: temp
    do i = 1, u%n
      y(i) = 0
    end do
    do j = 1, u%n
      temp = x(j)
      do i = 1, j
        y(i) = y(i) + temp*u%at(i, j)
      end do
    end do
  end subroutine

  pure real(dp) function packed_at(this, i, j)
    class(packed_upper), intent(in) :: this
    integer, intent(in) :: i, j
    packed_at = this%ap(i + (j*(j - 1))/2)
  end function

end module blas2_modern
