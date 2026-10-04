! CSR SpMV abstraction-penalty driver.
!
! Times the F77 baseline (spmv_f77.f) and each modern variant
! (spmv_modern.f90) on four matrices that span row length and working-set
! size.  Kernels are compiled separately and called without LTO, so nothing
! is inlined into the timing loops.
!
! Every variant reads the SAME storage: the matrix lives only in the
! derived-type object A, and the array-argument variants are passed its
! components.  The variant order is rotated every trial.  (With a separate
! copy for the derived-type variants and a fixed order, LLC-sized matrices
! showed up to 1.8x "penalties" for byte-identical code: whichever variant
! first switched copies ran cold.)
!
! Output, one line per (matrix, variant):
!   matrix  variant  n  nnz  t_f77  t_variant      [ns per call]
program bench_spmv
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use spmv_modern
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: nvar = 10, ntrial = 9, nmat = 4
  real(dp), parameter :: tmin_trial = 0.02_dp      ! seconds per trial
  character(len=10), parameter :: vname(nvar) = [character(len=10) :: &
       'f77', 'ctrl', 'ashape', 'dotprod', 'sum', 'assoc', 'dtype', 'tbp', 'rowfn', 'dc']
  character(len=10), parameter :: mname(nmat) = [character(len=10) :: &
       'lap2d_32', 'lap2d_1000', 'rand_64', 'rand_500']

  integer, allocatable :: rowptr(:), col(:)
  real(dp), allocatable :: val(:), x(:), y(:), yref(:)
  type(csr_matrix) :: a
  integer :: n, im, v, vv, reps, trial
  real(dp) :: best(nvar), t

  do im = 1, nmat
    select case (im)
    case (1); call make_lap2d(32)
    case (2); call make_lap2d(1000)
    case (3); call make_random(20000, 64)
    case (4); call make_random(2000, 500)
    end select
    call finish_setup()
    call check_equivalence(im)

    reps = 1
    do
      t = run(1, reps)
      if (t >= tmin_trial) exit
      reps = 2*reps
    end do
    best = huge(1.0_dp)
    do trial = 1, ntrial
      do vv = 0, nvar - 1
        v = 1 + mod(vv + trial, nvar)
        best(v) = min(best(v), run(v, reps))
      end do
    end do
    do v = 2, nvar
      write(*, '(a10,1x,a10,1x,i8,1x,i9,2(1x,es12.5))') mname(im), vname(v), n, &
           a%rowptr(n+1) - 1, 1.0e9_dp*best(1)/reps, 1.0e9_dp*best(v)/reps
    end do
  end do

contains

  !> 5-point Laplacian pattern on an m x m grid, random values.
  !> Rows have 3-5 entries: the short-row case.
  subroutine make_lap2d(m)
    integer, intent(in) :: m
    integer :: ix, iy, i, k
    n = m*m
    if (allocated(rowptr)) deallocate(rowptr, col, val)
    allocate(rowptr(n+1), col(5*n), val(5*n))
    k = 0
    rowptr(1) = 1
    do iy = 1, m
      do ix = 1, m
        i = ix + (iy - 1)*m
        if (iy > 1) then
          k = k + 1; col(k) = i - m
        end if
        if (ix > 1) then
          k = k + 1; col(k) = i - 1
        end if
        k = k + 1; col(k) = i
        if (ix < m) then
          k = k + 1; col(k) = i + 1
        end if
        if (iy < m) then
          k = k + 1; col(k) = i + m
        end if
        rowptr(i+1) = k + 1
      end do
    end do
    call random_number(val)
  end subroutine

  !> nrow x nrow matrix with exactly len entries per row, columns sorted and
  !> spread over the whole range (jittered regular spacing).
  subroutine make_random(nrow, len)
    integer, intent(in) :: nrow, len
    integer :: i, k, kk
    real(dp) :: u(len)
    n = nrow
    if (allocated(rowptr)) deallocate(rowptr, col, val)
    allocate(rowptr(n+1), col(n*len), val(n*len))
    kk = 0
    do i = 1, n
      rowptr(i) = kk + 1
      call random_number(u)
      do k = 1, len
        kk = kk + 1
        col(kk) = 1 + int((k - 1 + u(k))*real(n, dp)/len)
      end do
    end do
    rowptr(n+1) = kk + 1
    call random_number(val)
    val = 2*val - 1
  end subroutine

  subroutine finish_setup()
    if (allocated(x)) deallocate(x, y, yref)
    allocate(x(n), y(n), yref(n))
    call random_number(x)
    x = 2*x - 1
    ! single copy of the matrix, owned by A; work arrays released
    a%n = n
    a%col = col(1:rowptr(n+1)-1)
    a%val = val(1:rowptr(n+1)-1)
    call move_alloc(rowptr, a%rowptr)
    deallocate(col, val)
  end subroutine

  real(dp) function run(v, reps) result(elapsed)
    integer, intent(in) :: v, reps
    integer(int64) :: c0, c1, rate
    integer :: r
    call system_clock(c0, rate)
    select case (v)
    case (1); do r = 1, reps; call fcsrmv(n, a%rowptr, a%col, a%val, x, y); end do
    case (2); do r = 1, reps; call csrmv_ctrl(n, a%rowptr, a%col, a%val, x, y); end do
    case (3); do r = 1, reps; call csrmv_ashape(a%rowptr, a%col, a%val, x, y); end do
    case (4); do r = 1, reps; call csrmv_dotprod(n, a%rowptr, a%col, a%val, x, y); end do
    case (5); do r = 1, reps; call csrmv_sum(n, a%rowptr, a%col, a%val, x, y); end do
    case (6); do r = 1, reps; call csrmv_assoc(n, a%rowptr, a%col, a%val, x, y); end do
    case (7); do r = 1, reps; call csrmv_dtype(a, x, y); end do
    case (8); do r = 1, reps; call a%matvec(x, y); end do
    case (9); do r = 1, reps; call csrmv_rowfn(a, x, y); end do
    case (10); do r = 1, reps; call csrmv_dc(n, a%rowptr, a%col, a%val, x, y); end do
    case default; error stop 'unknown variant'
    end select
    call system_clock(c1)
    elapsed = real(c1 - c0, dp)/real(rate, dp)
  end function

  !> All variants do identical arithmetic, so results must match bitwise.
  subroutine check_equivalence(im)
    integer, intent(in) :: im
    integer :: v
    real(dp) :: t, rel
    t = run(1, 1)
    yref = y
    do v = 2, nvar
      y = 0
      t = run(v, 1)
      rel = maxval(abs(y - yref))/maxval(abs(yref))
      if (rel > 1.0e-12_dp) then
        write(*, '(a,a,1x,a,es10.2)') '# FAIL ', mname(im), vname(v), rel
        error stop 'variant does not reproduce the F77 result'
      else if (rel > 0) then
        write(*, '(a,a,1x,a,es10.2)') '# note: not bitwise: ', mname(im), vname(v), rel
      end if
    end do
  end subroutine

end program bench_spmv
