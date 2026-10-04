C     CSR sparse matrix-vector product Y := A*X, FORTRAN 77 style.
C     Row I holds entries VAL(K), column COLIDX(K), for
C     K = ROWPTR(I), ..., ROWPTR(I+1)-1.  Scalar accumulator, summation
C     in storage order.
      SUBROUTINE FCSRMV(N, ROWPTR, COLIDX, VAL, X, Y)
      INTEGER N, ROWPTR(N+1), COLIDX(*), I, K
      DOUBLE PRECISION VAL(*), X(*), Y(N), T
      DO 20 I = 1, N
         T = 0.0D0
         DO 10 K = ROWPTR(I), ROWPTR(I+1) - 1
            T = T + VAL(K)*X(COLIDX(K))
   10    CONTINUE
         Y(I) = T
   20 CONTINUE
      END
