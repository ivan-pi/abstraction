C     BLAS level-2 baselines, FORTRAN 77 style, after the reference
C     BLAS loop orderings (one UPLO/TRANS case each, ALPHA=1, BETA=0
C     unless stated).  Explicit-shape dummies with leading dimension LDA.
C
C     Y := A*X   (column / axpy form)
      SUBROUTINE FGEMVN(M, N, A, LDA, X, Y)
      INTEGER M, N, LDA, I, J
      DOUBLE PRECISION A(LDA,*), X(*), Y(*), TEMP
      DO 5 I = 1, M
         Y(I) = 0.0D0
    5 CONTINUE
      DO 20 J = 1, N
         TEMP = X(J)
         DO 10 I = 1, M
            Y(I) = Y(I) + TEMP*A(I,J)
   10    CONTINUE
   20 CONTINUE
      END
C
C     Y := A**T*X   (dot form)
      SUBROUTINE FGEMVT(M, N, A, LDA, X, Y)
      INTEGER M, N, LDA, I, J
      DOUBLE PRECISION A(LDA,*), X(*), Y(*), TEMP
      DO 20 J = 1, N
         TEMP = 0.0D0
         DO 10 I = 1, M
            TEMP = TEMP + A(I,J)*X(I)
   10    CONTINUE
         Y(J) = TEMP
   20 CONTINUE
      END
C
C     A := A + ALPHA*X*Y**T
      SUBROUTINE FGER(M, N, ALPHA, X, Y, A, LDA)
      INTEGER M, N, LDA, I, J
      DOUBLE PRECISION ALPHA, X(*), Y(*), A(LDA,*), TEMP
      DO 20 J = 1, N
         TEMP = ALPHA*Y(J)
         DO 10 I = 1, M
            A(I,J) = A(I,J) + X(I)*TEMP
   10    CONTINUE
   20 CONTINUE
      END
C
C     Solve U*X = B in place, U upper triangular, non-unit, column form.
      SUBROUTINE FTRSV(N, A, LDA, X)
      INTEGER N, LDA, I, J
      DOUBLE PRECISION A(LDA,*), X(*), TEMP
      DO 20 J = N, 1, -1
         X(J) = X(J)/A(J,J)
         TEMP = X(J)
         DO 10 I = J - 1, 1, -1
            X(I) = X(I) - TEMP*A(I,J)
   10    CONTINUE
   20 CONTINUE
      END
C
C     Y := A*X, A symmetric, upper triangle referenced.  One fused
C     pass per column: axpy into Y(1:J-1) and dot with X(1:J-1).
      SUBROUTINE FSYMV(N, A, LDA, X, Y)
      INTEGER N, LDA, I, J
      DOUBLE PRECISION A(LDA,*), X(*), Y(*), TEMP1, TEMP2
      DO 5 I = 1, N
         Y(I) = 0.0D0
    5 CONTINUE
      DO 20 J = 1, N
         TEMP1 = X(J)
         TEMP2 = 0.0D0
         DO 10 I = 1, J - 1
            Y(I) = Y(I) + TEMP1*A(I,J)
            TEMP2 = TEMP2 + A(I,J)*X(I)
   10    CONTINUE
         Y(J) = Y(J) + TEMP1*A(J,J) + TEMP2
   20 CONTINUE
      END
C
C     Y := U*X, U upper triangular in column-packed storage AP.
      SUBROUTINE FTPMV(N, AP, X, Y)
      INTEGER N, I, J, K, KK
      DOUBLE PRECISION AP(*), X(*), Y(*), TEMP
      DO 5 I = 1, N
         Y(I) = 0.0D0
    5 CONTINUE
      KK = 1
      DO 20 J = 1, N
         TEMP = X(J)
         K = KK
         DO 10 I = 1, J
            Y(I) = Y(I) + TEMP*AP(K)
            K = K + 1
   10    CONTINUE
         KK = KK + J
   20 CONTINUE
      END
