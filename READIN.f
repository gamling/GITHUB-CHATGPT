      SUBROUTINE READIN(sefp,sefm,sebp,sebm,seup,seum,omega,n)
      IMPLICIT NONE
      INTEGER i,j,n,n0,kounter
      PARAMETER(n0=1000)
      DOUBLE PRECISION dummy
      DOUBLE PRECISION sefp(n),sefm(n),sebp(n),sebm(n)
      DOUBLE PRECISION seup(n),seum(n)
      DOUBLE PRECISION omega(n)
      DOUBLE PRECISION omega0(n0), aux(n0)

      OPEN(UNIT=20,FILE='SEFP.dat',STATUS='UNKNOWN')
      kounter=0
      DO i=1,n0
      READ(20,*,END=95) omega0(i),aux(i)
      kounter=kounter+1
      END DO
 95   CONTINUE
      CLOSE(20)
      CALL GINDERPOL(omega0,kounter,omega,n,aux,sefp)

      OPEN(UNIT=20,FILE='SEFM.dat',STATUS='UNKNOWN')
      DO i=1,kounter
      READ(20,*,END=96) dummy,aux(i)
      END DO
 96   CONTINUE
      CLOSE(20)
      CALL GINDERPOL(omega0,kounter,omega,n,aux,sefm)

      OPEN(UNIT=20,FILE='SEBP.dat',STATUS='UNKNOWN')
      DO i=1,kounter
      READ(20,*,END=97) omega0(i),aux(i)
      END DO
 97   CONTINUE
      CLOSE(20)
      CALL GINDERPOL(omega0,kounter,omega,n,aux,sebp)

      OPEN(UNIT=20,FILE='SEBM.dat',STATUS='UNKNOWN')
      DO i=1,kounter
      READ(20,*,END=98) dummy,aux(i)
      END DO
 98   CONTINUE
      CLOSE(20)
      CALL GINDERPOL(omega0,kounter,omega,n,aux,sebm)
      
      OPEN(UNIT=20,FILE='SEUP.dat',STATUS='UNKNOWN')
      DO i=1,kounter
      READ(20,*,END=88) dummy,aux(i)
      END DO
 88   CONTINUE
      CLOSE(20)
      CALL GINDERPOL(omega0,kounter,omega,n,aux,seup)
      OPEN(UNIT=20,FILE='SEUM.dat',STATUS='UNKNOWN')
      DO i=1,kounter
      READ(20,*,END=78) dummy,aux(i)
      END DO
 78   CONTINUE
      CLOSE(20)
      CALL GINDERPOL(omega0,kounter,omega,n,aux,seum)

      RETURN
      END
***********************************************************************
***********************************************************************
      SUBROUTINE wrpropa(fname,n,freq,repa,impa)
      IMPLICIT NONE
      INTEGER i,n
      DOUBLE PRECISION freq(1:n),repa(1:n),impa(1:n)
      CHARACTER(len=*) fname

      OPEN(UNIT=88,FILE=fname,STATUS='unknown')
      DO i=1,n
         WRITE(88,*) freq(i),repa(i),impa(i)
      END DO
      CLOSE(88)
      END
***********************************************************************
***********************************************************************
      SUBROUTINE wrout(fname,n,freq)
      IMPLICIT NONE
      INTEGER i,n
      DOUBLE PRECISION freq(1:n)
      CHARACTER(len=*) fname

      OPEN(UNIT=88,FILE=fname,STATUS='unknown')
      DO i=1,n
         WRITE(88,*) i,freq(i)
      END DO
      CLOSE(88)
      END
***********************************************************************
***********************************************************************
      SUBROUTINE rout(fname,n,freq,field)
      IMPLICIT NONE
      INTEGER i,n
      DOUBLE PRECISION freq(1:n),field(1:n)
      CHARACTER(len=*) fname

      OPEN(UNIT=88,FILE=fname,STATUS='unknown')
      DO i=1,n
         WRITE(88,*) freq(i),field(i)
      END DO
      CLOSE(88)
      END
***********************************************************************
***********************************************************************
