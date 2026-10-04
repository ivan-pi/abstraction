! BLAS level-2 abstraction-penalty driver.
!
! For each kernel, times the F77 baseline (blas2_f77.f) and each modern
! variant (blas2_modern.f90, blas2_pdt.f90) on square n x n problems stored
! with leading dimension lda = n + 4, so column sections are not contiguous.
! Kernels are compiled separately and called without LTO, so nothing is
! inlined into the timing loops.
!
! Output, one line per (kernel, variant, n):
!   kernel  variant  n  t_f77  t_variant      [ns per call]
program bench2
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use blas2_modern
#ifdef HAVE_PDT
  use blas2_pdt
#endif
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: nk = 6, maxv = 7, ntrial = 9, pad = 4
  integer, parameter :: sizes(4) = [8, 32, 256, 2000]
  real(dp), parameter :: tmin_trial = 0.02_dp      ! seconds per trial
  real(dp) :: alpha = 1.0e-12_dp     ! tiny while timing (keeps repeated GER bounded)

  character(len=6), parameter :: kname(nk) = [character(len=6) :: &
       'gemv_n', 'gemv_t', 'ger', 'trsv', 'symv', 'tpmv']
  character(len=10) :: vname(maxv, nk)
  integer :: nvar(nk)

  real(dp), allocatable :: a(:,:), ap(:), x(:), b(:), y(:), yref(:), aref(:,:), gref(:,:)
  type(dmat) :: dm
  type(packed_upper) :: pu
#ifdef HAVE_PDT
  type(pmat(:,:)), allocatable :: pm
#endif
  integer :: n, lda, is, k, v, reps, trial
  real(dp) :: best(maxv), t

  call setup_names()

  do is = 1, size(sizes)
    n = sizes(is)
    lda = n + pad
    call init(n, lda)
    call check_equivalence(n)
    do k = 1, nk
      reps = 1
      do
        t = run(k, 1, reps)
        if (t >= tmin_trial) exit
        reps = 2*reps
      end do
      best = huge(1.0_dp)
      do trial = 1, ntrial
        do v = 1, nvar(k)
          best(v) = min(best(v), run(k, v, reps))
        end do
      end do
      do v = 2, nvar(k)
        write(*, '(a6,1x,a10,1x,i5,2(1x,es12.5))') kname(k), vname(v, k), n, &
             1.0e9_dp*best(1)/reps, 1.0e9_dp*best(v)/reps
      end do
    end do
  end do

contains

  subroutine setup_names()
    vname = ''
    vname(1:7, 1) = [character(len=10) :: 'f77', 'arrsec', 'ashape', 'matmul', &
                     'dtype', 'tbp', 'pdt']
    nvar(1) = 6
#ifdef HAVE_PDT
    nvar(1) = 7
#endif
    vname(1:5, 2) = [character(len=10) :: 'f77', 'sum', 'dotprod', 'matmul_vm', 'matmul_tr']
    nvar(2) = 5
    vname(1:5, 3) = [character(len=10) :: 'f77', 'arrsec', 'dc_ji', 'dc_ij', 'spread']
    nvar(3) = 5
    vname(1:3, 4) = [character(len=10) :: 'f77', 'arrsec', 'arrsec_tmp']
    nvar(4) = 3
    vname(1:3, 5) = [character(len=10) :: 'f77', 'arrsec_sum', 'dotprod']
    nvar(5) = 3
    vname(1:3, 6) = [character(len=10) :: 'f77', 'arrsec', 'accessor']
    nvar(6) = 3
  end subroutine

  !> Random A in [-1,1] with a dominant diagonal (keeps TRSV tame), packed
  !> copy of its upper triangle, and the derived-type copies.
  subroutine init(n, lda)
    integer, intent(in) :: n, lda
    integer :: i, j, kk
    if (allocated(a)) deallocate(a, ap, x, b, y, yref, aref)
    allocate(a(lda, n), aref(lda, n), ap(n*(n+1)/2), x(n), b(n), y(n), yref(n))
    call random_number(a)
    a = 2*a - 1
    do j = 1, n
      a(j, j) = a(j, j) + n
    end do
    call random_number(b)
    b = 2*b - 1
    x = b
    kk = 0
    do j = 1, n
      do i = 1, j
        kk = kk + 1
        ap(kk) = a(i, j)
      end do
    end do
    dm%m = n
    dm%n = n
    dm%a = a
    pu%n = n
    pu%ap = ap
#ifdef HAVE_PDT
    if (allocated(pm)) deallocate(pm)
    allocate(pmat(lda, n) :: pm)
    pm%m = n
    pm%a = a
#endif
  end subroutine

  !> Time `reps` calls of kernel k, variant v.  TRSV re-copies the RHS
  !> before every solve (O(n) next to the O(n^2) solve, same for all).
  real(dp) function run(k, v, reps) result(elapsed)
    integer, intent(in) :: k, v, reps
    integer(int64) :: c0, c1, rate
    integer :: r

    call system_clock(c0, rate)
    select case (100*k + v)
    ! GEMV N
    case (101); do r = 1, reps; call fgemvn(n, n, a, lda, x, y); end do
    case (102); do r = 1, reps; call gemvn_arrsec(n, n, a, lda, x, y); end do
    case (103); do r = 1, reps; call gemvn_ashape(a(1:n,1:n), x, y); end do
    case (104); do r = 1, reps; call gemvn_matmul(n, n, a, lda, x, y); end do
    case (105); do r = 1, reps; call gemvn_dtype(dm, x, y); end do
    case (106); do r = 1, reps; call dm%mv(x, y); end do
#ifdef HAVE_PDT
    case (107); do r = 1, reps; call gemvn_pdt(pm, x, y); end do
#endif
    ! GEMV T
    case (201); do r = 1, reps; call fgemvt(n, n, a, lda, x, y); end do
    case (202); do r = 1, reps; call gemvt_sum(n, n, a, lda, x, y); end do
    case (203); do r = 1, reps; call gemvt_dotprod(n, n, a, lda, x, y); end do
    case (204); do r = 1, reps; call gemvt_matmul_vm(n, n, a, lda, x, y); end do
    case (205); do r = 1, reps; call gemvt_matmul_tr(n, n, a, lda, x, y); end do
    ! GER
    case (301); do r = 1, reps; call fger(n, n, alpha, x, b, a, lda); end do
    case (302); do r = 1, reps; call ger_arrsec(n, n, alpha, x, b, a, lda); end do
    case (303); do r = 1, reps; call ger_dc_ji(n, n, alpha, x, b, a, lda); end do
    case (304); do r = 1, reps; call ger_dc_ij(n, n, alpha, x, b, a, lda); end do
    case (305); do r = 1, reps; call ger_spread(n, n, alpha, x, b, a, lda); end do
    ! TRSV
    case (401); do r = 1, reps; y = b; call ftrsv(n, a, lda, y); end do
    case (402); do r = 1, reps; y = b; call trsv_arrsec(n, a, lda, y); end do
    case (403); do r = 1, reps; y = b; call trsv_arrsec_tmp(n, a, lda, y); end do
    ! SYMV
    case (501); do r = 1, reps; call fsymv(n, a, lda, x, y); end do
    case (502); do r = 1, reps; call symv_arrsec_sum(n, a, lda, x, y); end do
    case (503); do r = 1, reps; call symv_dotprod(n, a, lda, x, y); end do
    ! TPMV
    case (601); do r = 1, reps; call ftpmv(n, ap, x, y); end do
    case (602); do r = 1, reps; call tpmv_arrsec(n, ap, x, y); end do
    case (603); do r = 1, reps; call tpmv_accessor(pu, x, y); end do
    case default; error stop 'unknown kernel/variant'
    end select
    call system_clock(c1)
    elapsed = real(c1 - c0, dp)/real(rate, dp)
  end function

  !> Every variant must reproduce the F77 result.  Report bitwise equality
  !> or the relative difference; fail beyond 1e-12.  GER is checked with
  !> alpha = 0.5 so the update is not lost in rounding.
  subroutine check_equivalence(n)
    integer, intent(in) :: n
    integer :: k, v
    real(dp) :: rel, t
    alpha = 0.5_dp
    aref = a
    do k = 1, nk
      do v = 1, nvar(k)
        t = run(k, v, 1)
        if (v == 1) then
          call save_ref(k)
        else
          rel = compare(k)
          if (rel > 1.0e-12_dp) then
            write(*, '(a,a,1x,a,i6,es10.2)') '# FAIL ', kname(k), vname(v, k), n, rel
            error stop 'variant does not reproduce the F77 result'
          else if (rel > 0) then
            write(*, '(a,a,1x,a,i6,es10.2)') '# note: not bitwise: ', kname(k), &
                 vname(v, k), n, rel
          end if
        end if
        a = aref
      end do
    end do
    alpha = 1.0e-12_dp
  end subroutine

  ! The reference result: y for every kernel except GER, whose result is
  ! the updated matrix.
  subroutine save_ref(k)
    integer, intent(in) :: k
    if (k == 3) then
      gref = a(1:n,1:n)
    else
      yref = y
    end if
  end subroutine

  real(dp) function compare(k) result(rel)
    integer, intent(in) :: k
    if (k == 3) then
      rel = maxval(abs(a(1:n,1:n) - gref))/maxval(abs(gref))
    else
      rel = maxval(abs(y - yref))/maxval(abs(yref))
    end if
  end function

end program bench2
