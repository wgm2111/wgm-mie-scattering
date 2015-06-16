C     file : matr_.f
C 
c$$$
c$$$      PROGRAM MATR_
c$$$C     
c$$$C     Purpose: 
c$$$C       Test the routine MATR before wrapping with f2py 
c$$$C       for use in python.
c$$$C
c$$$      
c$$$      PARAMETER (NANG=181, NAB=600)
c$$$     
c$$$      IMPLICIT REAL*8 (A-H,O-Z)                                                 
c$$$
c$$$      INTEGER i
c$$$      REAL*8 ANG(NANG), F11(NANG), F33(NANG), F12(NANG), F34(NANG) 
c$$$      REAL*8 A1(NAB), A2(NAB), A3(NAB), A4(NAB), B1(NAB), B2(NAB)
c$$$
c$$$
c$$$      PI = DACOS(-1D0)
c$$$
c$$$c     Initialize the input variables 
c$$$      DO i = 1, NANG
c$$$         ANG(i) = PI * (DFLOAT(i)-1D0) / (NANG-1)
c$$$      END DO
c$$$
c$$$      F11 = 0.0
c$$$      F33 = 0.0
c$$$      F12 = 0.0
c$$$      F34 = 0.0
c$$$
c$$$c     Initialize the output variables
c$$$      
c$$$      A1 = (/ (1D0 / DFLOAT(i), i=1,NAB) /)
c$$$      A2 = (/ (1D0 / DFLOAT(i)**2, i=1,NAB) /)
c$$$      A3 = (/ (1D0 / DFLOAT(i)**2, i=1,NAB) /)
c$$$      A4 = (/ (1D0 / DFLOAT(i)**2, i=1,NAB) /)
c$$$      B1 = (/ (1D0 / DFLOAT(i)**2, i=1,NAB) /)
c$$$      B2 = (/ (1D0 / DFLOAT(i)**2, i=1,NAB) /)
c$$$
c$$$c     Call the matr subroutine to compute scattering matrix
c$$$      CALL MATR(ANG,NANG, A1,A2,A3,A4,B1,B2,NAB, F11,F33,F12,F34)
c$$$
c$$$c     Print the input values
c$$$c      WRITE(*,2001) 'A1', 'A2', 'A3', 'A4', 'B1', 'B2'
c$$$c 2001 FORMAT(1X, 6A12)      
c$$$c      DO i = 1, NAB
c$$$c         WRITE(*,2002) A1(i), A2(i), A3(i), A4(i), B1(i), B2(i)
c$$$c      END DO
c$$$c 2002 FORMAT(1X, 6E12.4)
c$$$
c$$$c     Print the output variables 
c$$$c      WRITE(*,1001) 'ANGLE', 'F11', 'F33', 'F12', 'F34'
c$$$c 1001 FORMAT(1X, 5A16)      
c$$$c      DO i = 1, NANG
c$$$c         WRITE(*,1002) ANG(i), F11(i), F33(i), F12(i), F34(i) 
c$$$c      END DO
c$$$c 1002 FORMAT(1X, 5E16.4)
c$$$
c$$$      END PROGRAM MATR_



C****************************************************************
 
C     CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
C     COEFFICIENTS
 
C     A1,...,B2 - EXPANSION COEFFICIENTS
C     NAB       - number of expansio coefficients


C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
C    NPNA - NUMBER OF SCATTERING ANGLES
C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
C        180*(I-1)/(NPNA-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
      SUBROUTINE MATR (ANG,NANG, A1,A2,A3,A4,B1,B2,NAB, F11array, 
     &     F33array, F12array, F34array)

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NANG, NAB
      REAL*8 ANG(NANG), F11array(NANG), F33array(NANG), F12array(NANG),
     &     F34array(NANG)
      REAL*8 A1(NAB),A2(NAB),A3(NAB),A4(NAB),B1(NAB),B2(NAB)
Cf2py intent(in) ANG, A1, A2, A3, A4, B1, B2
Cf2py intent(out) F11array, F33array, F12array, F34array
      LMAX = NAB-1
      L1MAX = NAB
      NPNA = NANG
      D6=DSQRT(6D0)*0.25D0
      N=NPNA
      PI = DACOS(-1D0)

      DO 500 I1=1,N
         TAA= ANG(I1)         
         TB= ANG(I1) * PI / 180D0
         U=DCOS(TAA)
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1)*DFLOAT(L*L-4)
            PL4=1D0/(DFLOAT(L)*DFLOAT(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400    CONTINUE

         F33=(F2-F3)*0.5D0
         
C        Store current values of scattering matrix for output
         F11array(I1) = F11
         F33array(I1) = F33
         F12array(I1) = F12
         F34array(I1) = F34
 500  CONTINUE
      RETURN
      END
 
