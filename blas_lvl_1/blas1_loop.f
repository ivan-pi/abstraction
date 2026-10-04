C     BLAS level-1 kernels, FORTRAN 77 style: explicit-shape dummies,
C     labelled DO loops, scalar temporaries.  Unit stride unless the
C     routine takes INCX/INCY.  No reassociation-sensitive tricks
C     (no unrolling, no partial sums): the compiler sees the plain loop.
C
      SUBROUTINE LCOPY(N, X, Y)
      INTEGER N, I
      DOUBLE PRECISION X(N), Y(N)
      DO 10 I = 1, N
         Y(I) = X(I)
   10 CONTINUE
      END
C
      SUBROUTINE LSCAL(N, A, X)
      INTEGER N, I
      DOUBLE PRECISION A, X(N)
      DO 10 I = 1, N
         X(I) = A*X(I)
   10 CONTINUE
      END
C
      SUBROUTINE LAXPY(N, A, X, Y)
      INTEGER N, I
      DOUBLE PRECISION A, X(N), Y(N)
      DO 10 I = 1, N
         Y(I) = Y(I) + A*X(I)
   10 CONTINUE
      END
C
      SUBROUTINE LAXPBY(N, A, X, B, Y)
      INTEGER N, I
      DOUBLE PRECISION A, B, X(N), Y(N)
      DO 10 I = 1, N
         Y(I) = A*X(I) + B*Y(I)
   10 CONTINUE
      END
C
      SUBROUTINE LWAXPB(N, A, X, B, Y, W)
      INTEGER N, I
      DOUBLE PRECISION A, B, X(N), Y(N), W(N)
      DO 10 I = 1, N
         W(I) = A*X(I) + B*Y(I)
   10 CONTINUE
      END
C
      DOUBLE PRECISION FUNCTION LDOT(N, X, Y)
      INTEGER N, I
      DOUBLE PRECISION X(N), Y(N), S
      S = 0.0D0
      DO 10 I = 1, N
         S = S + X(I)*Y(I)
   10 CONTINUE
      LDOT = S
      END
C
      DOUBLE PRECISION FUNCTION LASUM(N, X)
      INTEGER N, I
      DOUBLE PRECISION X(N), S
      S = 0.0D0
      DO 10 I = 1, N
         S = S + ABS(X(I))
   10 CONTINUE
      LASUM = S
      END
C
C     Naive (unscaled) 2-norm, same algorithm as SQRT(SUM(X**2)).
      DOUBLE PRECISION FUNCTION LNRM2(N, X)
      INTEGER N, I
      DOUBLE PRECISION X(N), S
      S = 0.0D0
      DO 10 I = 1, N
         S = S + X(I)*X(I)
   10 CONTINUE
      LNRM2 = SQRT(S)
      END
C
C     First index of max |x(i)|, same tie rule as MAXLOC.
      INTEGER FUNCTION LIAMAX(N, X)
      INTEGER N, I, K
      DOUBLE PRECISION X(N), XMAX
      K = 1
      XMAX = ABS(X(1))
      DO 10 I = 2, N
         IF (ABS(X(I)) .GT. XMAX) THEN
            K = I
            XMAX = ABS(X(I))
         END IF
   10 CONTINUE
      LIAMAX = K
      END
C
C     Strided AXPY, reference-BLAS indexing (positive increments).
      SUBROUTINE LAXPYI(N, A, X, INCX, Y, INCY)
      INTEGER N, INCX, INCY, I, IX, IY
      DOUBLE PRECISION A, X(*), Y(*)
      IX = 1
      IY = 1
      DO 10 I = 1, N
         Y(IY) = Y(IY) + A*X(IX)
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      END
C
      DOUBLE PRECISION FUNCTION LDOTI(N, X, INCX, Y, INCY)
      INTEGER N, INCX, INCY, I, IX, IY
      DOUBLE PRECISION X(*), Y(*), S
      S = 0.0D0
      IX = 1
      IY = 1
      DO 10 I = 1, N
         S = S + X(IX)*Y(IY)
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      LDOTI = S
      END
