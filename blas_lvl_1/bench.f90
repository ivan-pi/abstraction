! Driver: times each BLAS-1 kernel in loop form (L*) and array form (A*).
! Kernels live in separately compiled files and are called as external
! procedures, so nothing is inlined into the timing loop (build without LTO).
!
! Output, one line per (kernel, n):  name  n  t_loop  t_array   [ns per call]
program bench
  use, intrinsic :: iso_fortran_env, only: int64, real64
  implicit none

  integer, parameter :: nk = 11, ntrial = 11
  character(len=8), parameter :: names(nk) = [character(len=8) :: &
       'copy', 'scal', 'axpy', 'axpby', 'waxpby', 'dot', 'asum', 'nrm2', &
       'iamax', 'axpy_s2', 'dot_s2']
  integer, parameter :: sizes(3) = [1000, 50000, 5000000]
  real(real64), parameter :: tmin_trial = 0.02_real64   ! seconds per trial

  double precision, external :: ldot, lasum, lnrm2, ldoti
  double precision, external :: adot, aasum, anrm2, adoti
  integer, external :: liamax, aiamax

  real(real64), allocatable :: x(:), y(:), w(:)
  real(real64) :: sink, tl, ta, best(2), t
  integer :: is, k, n, reps, trial, v

  allocate(x(2*maxval(sizes)), y(2*maxval(sizes)), w(2*maxval(sizes)))
  sink = 0

  call check_equivalence()

  do is = 1, size(sizes)
    n = sizes(is)
    do k = 1, nk
      call init()
      ! calibrate repetitions on the loop variant
      reps = 1
      do
        t = run(k, 1, n, reps)
        if (t >= tmin_trial) exit
        reps = reps*2
      end do
      ! alternate variants to cancel drift; keep the minimum
      best = huge(1.0_real64)
      do trial = 1, ntrial
        do v = 1, 2
          best(v) = min(best(v), run(k, v, n, reps))
        end do
      end do
      tl = 1.0e9_real64*best(1)/reps
      ta = 1.0e9_real64*best(2)/reps
      write(*, '(a8,1x,i9,2(1x,es12.5))') names(k), n, tl, ta
    end do
  end do
  write(*, '(a,es12.4)') '# sink ', sink

contains

  subroutine init()
    call random_number(x)
    call random_number(y)
    x = 2*x - 1
    y = 2*y - 1
    w = 0
  end subroutine

  ! Time `reps` calls of kernel k in variant v (1 = loop, 2 = array).
  real(real64) function run(k, v, n, reps) result(elapsed)
    integer, intent(in) :: k, v, n, reps
    integer(int64) :: c0, c1, rate
    integer :: r
    real(real64), parameter :: a = 1.0e-12_real64, b = 1.0_real64 - 1.0e-12_real64
    real(real64), parameter :: s = 1.0_real64 + 1.0e-12_real64

    call system_clock(c0, rate)
    select case (k*10 + v)
    case (11); do r = 1, reps; call lcopy(n, x, y); end do
    case (12); do r = 1, reps; call acopy(n, x, y); end do
    case (21); do r = 1, reps; call lscal(n, s, x); end do
    case (22); do r = 1, reps; call ascal(n, s, x); end do
    case (31); do r = 1, reps; call laxpy(n, a, x, y); end do
    case (32); do r = 1, reps; call aaxpy(n, a, x, y); end do
    case (41); do r = 1, reps; call laxpby(n, a, x, b, y); end do
    case (42); do r = 1, reps; call aaxpby(n, a, x, b, y); end do
    case (51); do r = 1, reps; call lwaxpb(n, a, x, b, y, w); end do
    case (52); do r = 1, reps; call awaxpb(n, a, x, b, y, w); end do
    case (61); do r = 1, reps; sink = sink + ldot(n, x, y); end do
    case (62); do r = 1, reps; sink = sink + adot(n, x, y); end do
    case (71); do r = 1, reps; sink = sink + lasum(n, x); end do
    case (72); do r = 1, reps; sink = sink + aasum(n, x); end do
    case (81); do r = 1, reps; sink = sink + lnrm2(n, x); end do
    case (82); do r = 1, reps; sink = sink + anrm2(n, x); end do
    case (91); do r = 1, reps; sink = sink + liamax(n, x); end do
    case (92); do r = 1, reps; sink = sink + aiamax(n, x); end do
    case (101); do r = 1, reps; call laxpyi(n, a, x, 2, y, 2); end do
    case (102); do r = 1, reps; call aaxpyi(n, a, x, 2, y, 2); end do
    case (111); do r = 1, reps; sink = sink + ldoti(n, x, 2, y, 2); end do
    case (112); do r = 1, reps; sink = sink + adoti(n, x, 2, y, 2); end do
    end select
    call system_clock(c1)
    elapsed = real(c1 - c0, real64)/real(rate, real64)
  end function

  ! Both variants must compute the same thing.  Elementwise kernels and
  ! iamax must match exactly; reductions may differ only by summation order
  ! (relevant under -ffast-math or similar reassociation).
  subroutine check_equivalence()
    integer, parameter :: m = 1003
    real(real64) :: x1(2*m), y1(2*m), y2(2*m), w1(m), w2(m), r1, r2, err

    call random_number(x1); x1 = 2*x1 - 1
    call random_number(y1); y1 = 2*y1 - 1
    err = 0

    y2 = y1; call lcopy(m, x1, y2); w1 = y2(1:m)
    y2 = y1; call acopy(m, x1, y2); err = max(err, maxval(abs(w1 - y2(1:m))))
    y2 = y1; call lscal(m, 1.5d0, y2); w1 = y2(1:m)
    y2 = y1; call ascal(m, 1.5d0, y2); err = max(err, maxval(abs(w1 - y2(1:m))))
    y2 = y1; call laxpy(m, 0.7d0, x1, y2); w1 = y2(1:m)
    y2 = y1; call aaxpy(m, 0.7d0, x1, y2); err = max(err, maxval(abs(w1 - y2(1:m))))
    y2 = y1; call laxpby(m, 0.7d0, x1, 0.3d0, y2); w1 = y2(1:m)
    y2 = y1; call aaxpby(m, 0.7d0, x1, 0.3d0, y2); err = max(err, maxval(abs(w1 - y2(1:m))))
    call lwaxpb(m, 0.7d0, x1, 0.3d0, y1, w1)
    call awaxpb(m, 0.7d0, x1, 0.3d0, y1, w2); err = max(err, maxval(abs(w1 - w2)))
    y2 = y1; call laxpyi(m, 0.7d0, x1, 2, y2, 2); w1 = y2(1:2*m:2)
    y2 = y1; call aaxpyi(m, 0.7d0, x1, 2, y2, 2); err = max(err, maxval(abs(w1 - y2(1:2*m:2))))
    if (err /= 0) error stop 'elementwise kernels differ'

    r1 = ldot(m, x1, y1);  r2 = adot(m, x1, y1);  call chk('dot', r1, r2)
    r1 = lasum(m, x1);     r2 = aasum(m, x1);     call chk('asum', r1, r2)
    r1 = lnrm2(m, x1);     r2 = anrm2(m, x1);     call chk('nrm2', r1, r2)
    r1 = ldoti(m, x1, 2, y1, 2); r2 = adoti(m, x1, 2, y1, 2); call chk('dot_s2', r1, r2)
    if (liamax(m, x1) /= aiamax(m, x1)) error stop 'iamax differs'
    write(*, '(a)') '# equivalence check passed'
  end subroutine

  subroutine chk(name, p, q)
    character(*), intent(in) :: name
    real(real64), intent(in) :: p, q
    real(real64), parameter :: tol = 1.0e-12_real64
    if (abs(p - q) > tol*max(1.0_real64, abs(p))) then
      write(*, *) name, p, q
      error stop 'reduction kernels differ'
    end if
    if (p /= q) write(*, '(3a)') '# note: ', name, ' differs in last bits (summation order)'
  end subroutine

end program
