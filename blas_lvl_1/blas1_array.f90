! BLAS level-1 kernels, array-expression style.
! Interfaces are identical to blas1_loop.f (explicit-shape / assumed-size
! dummies, external procedures), so the only difference the compiler sees
! is DO loop vs array expression / transformational intrinsic.

subroutine acopy(n, x, y)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: y(n)
  y = x
end subroutine

subroutine ascal(n, a, x)
  integer, intent(in) :: n
  double precision, intent(in) :: a
  double precision, intent(inout) :: x(n)
  x = a*x
end subroutine

subroutine aaxpy(n, a, x, y)
  integer, intent(in) :: n
  double precision, intent(in) :: a, x(n)
  double precision, intent(inout) :: y(n)
  y = y + a*x
end subroutine

subroutine aaxpby(n, a, x, b, y)
  integer, intent(in) :: n
  double precision, intent(in) :: a, b, x(n)
  double precision, intent(inout) :: y(n)
  y = a*x + b*y
end subroutine

subroutine awaxpb(n, a, x, b, y, w)
  integer, intent(in) :: n
  double precision, intent(in) :: a, b, x(n), y(n)
  double precision, intent(out) :: w(n)
  w = a*x + b*y
end subroutine

double precision function adot(n, x, y)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n), y(n)
  adot = sum(x*y)
end function

double precision function aasum(n, x)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  aasum = sum(abs(x))
end function

double precision function anrm2(n, x)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  anrm2 = sqrt(sum(x*x))
end function

integer function aiamax(n, x)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  aiamax = maxloc(abs(x), dim=1)
end function

subroutine aaxpyi(n, a, x, incx, y, incy)
  integer, intent(in) :: n, incx, incy
  double precision, intent(in) :: a, x(*)
  double precision, intent(inout) :: y(*)
  y(1:1+(n-1)*incy:incy) = y(1:1+(n-1)*incy:incy) + a*x(1:1+(n-1)*incx:incx)
end subroutine

double precision function adoti(n, x, incx, y, incy)
  integer, intent(in) :: n, incx, incy
  double precision, intent(in) :: x(*), y(*)
  adoti = sum(x(1:1+(n-1)*incx:incx)*y(1:1+(n-1)*incy:incy))
end function
