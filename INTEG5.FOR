
      SUBROUTINE INTEG5(isym,tsum)
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical isym(2)
c
      parameter(eta1 = 1.0e-16, eta3 = 10.0)
c
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(4),ZA(4),FM(4),FMA(4),RC(4),AC(4),PNLOC(4),FS(4),ECM(4)
     2,FK(4),FK2(4),ETA(4),DR(4),FMU(4),ACC(8),DELQ,AMPL(8)
     3,K,KZ,LPLUS,LPL2,IS(4),NS(4),LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,LTRP(8),JTRP(8),ISTRP(8)
     5,MLTR,IBUFF,IWORD,KMAX,NMAX,LMAX,IFLL,ICO,is_off,iextra
c
      Common/array0/Space0(ispc0)
      Common/array1/Space1(ispc1)
      Common/array2/Space2(ispc2)
      Common/array3/Space3(ispc3)
      Common/n1n2ez/n1(nl2),n2(nl2),ezer(jmt)
c
      dimension u(1600),v(1600)
     1,         f(nl2),fp(nl2),g(nl2),gp(nl2),s(nl2)
     2,         amat(imat*imat),cmat(16*jmt),e(16*nl2),wait(1600)
     3,N3(nl2),N4(nl2)        ,UZ(jmt),R1(jmt),R2(jmt)
     4,B1(jmt),B2(jmt),P1(jmt),P2(jmt),S1(jmt),S2(jmt),S3(jmt),S4(jmt)
     5,C1(jmt),C2(jmt),C3(jmt),C4(jmt),D1(jmt),D2(jmt),D3(jmt),D4(jmt)
     6,aex(512)
c
      DIMENSION FACT(3),F1(6),F2(6),Q1(6),Q2(6),Q3(6),AKEX(2),I1MIN(2)
     1,I1MAX(2),WT(5),TEST(imat),AUG(2,2),A(2),B(2),Q(4),TSUM(2),T(2)
     2,CX(6),ANORM(6),SUM(6),AFK(2),ITEM(2)
     3,ETEMP(128),FF(2400)
c
c
      equivalence (Space0(   1),u(1)), (Space0(1601),v(1))
     1,           (Space0(3201      ),f(1)), (Space0(3201+  nl2),fp(1))
     2,           (Space0(3201+2*nl2),g(1)), (Space0(3201+3*nl2),gp(1))
     3,           (Space0(3201+4*nl2),s(1))
     4,           (Space0(3201+5*nl2),n3(1)),(Space0(3201+6*nl2),n4(1))
     4,           (Space1(   1),e(1)),    (Space1(3201),wait(1))
     5,           (Space1(4801),etemp(1)),(Space1(5001),aex(1))
      equivalence
     1       (Space1(5513       ),r1(1)),   (Space1(5513+ 1*jmt),r2(1))
     2,      (Space1(5513+ 2*jmt),b1(1)),   (Space1(5513+ 3*jmt),b2(1))
     3,      (Space1(5513+ 4*jmt),p1(1)),   (Space1(5513+ 5*jmt),p2(1))
     4,      (Space1(5513+ 6*jmt),s1(1)),   (Space1(5513+ 7*jmt),s2(1))
     5,      (Space1(5513+ 8*jmt),s3(1)),   (Space1(5513+ 9*jmt),s4(1))
     6,      (Space1(5513+10*jmt),c1(1)),   (Space1(5513+11*jmt),c2(1))
     7,      (Space1(5513+12*jmt),c3(1)),   (Space1(5513+13*jmt),c4(1))
     8,      (Space1(5513+14*jmt),d1(1)),   (Space1(5513+15*jmt),d2(1))
     9,      (Space1(5513+16*jmt),d3(1)),   (Space1(5513+17*jmt),d4(1))
      equivalence
     1            (Space2(   1),amat(1)), (Space2(imat*imat+1),cmat(1))
     2,           (Space2(   1),FF(1)), (AKEX(1),ACC(7))
c
      DATA PI/3.1415926535898/
      DATA ABSC, WIDTH, DETZ, SCALE/ 8.0, 5.0, 1.E-6, 1.0E34/
C
      FUN(Zz)=SQRT( MAX ((Zz-DFACT**2)*ACCU**2+DFACT**2,0.D0))
C
      WRITE(6,9999)ALPHA
C      CALL KSECOND(TIME)
C      WRITE(6,9950)TIME
C     ACC1 CONTROLS WAVE NO. SPREAD FOR WAVE 1
      IF(ACC(1).EQ.0.0) ACC(1)=1.0
C     ACC2 CONTROLS WAVE NO. SPREAD FOR WAVE 2
      IF(ACC(2).EQ.0.0) ACC(2)=1.0
C     ACC3 CONTROLS WIDTH PARAMETER, GAM FOR WAVE 1
C     ACC4 CONTROLS WIDTH PARAMETER, GAM FOR WAVE 2
C     ACC5 CONTROLS TURNING POINT RADIUS FOR WAVE 1
C     ACC6 CONTROLS TURNING POINT RADIUS FOR WAVE 2
C     ACC7 CONTROLS CENTRAL K FOR WAVE NO. 1
      IF(ACC(7).EQ.0.0) ACC(7)=1.0
C     ACC8 CONTROLS CENTRAL K FOR WAVE NO. 2
      IF(ACC(8).EQ.0.0) ACC(8)=1.0
C
      kt = k
      WT(1)=375./288.
      WT(4)=WT(1)
      WT(2)=250./288.
      WT(3)=WT(2)
      WT(5)=190./288.
      LPL2=LPLUS+LPLUS
      kt=5*(kt/5)
      IETA3=ETA3
      RMAX=DRF*FLOAT(kt)
      DR2=DRF**2/12.0
C
      KMX2=2*KMAX
      DO 10 M=1,kt
      index = mod(m-1,5)+1
      MK=M+M-1
      WAIT(MK       )=WT(INDEX)
      WAIT(MK+1     )=0.0
      WAIT(MK  +KMX2)=WT(INDEX)
      WAIT(MK+1+KMX2)=0.0
   10 CONTINUE
      WAIT(MK       )=WAIT(MK       )*0.5
      WAIT(MK  +KMX2)=WAIT(MK  +KMX2)*0.5
C
      CALL FNLOC5(U(     1),V(    1),WAIT(      1)
     1           ,PNLOC(1),FK2(1),FK(1),ETA(1),RC(1),DRF,kt)
      CALL FNLOC5(U(KMX2+1),V(KMX2+1),WAIT(KMX2+1)
     1           ,PNLOC(2),FK2(2),FK(2),ETA(2),RC(2),DRF,kt)
C
      IL=0
      IX=0
      IZ=0
      DO 145 N=1,2
      ITEM(N)=0
      ACCU=ACC(N)
      ACCL=ACC(N+4)
C
C     SELECT WAVE NO. SPACING
C
C     CALCULATE ASYMPTOTIC WAVE NUMBER
C
      DFK=PI/RMAX
      fk2Rmax=FK2(N)-2.0*ETA(N)*FK(N)/RMAX
      if(fk2Rmax.gt.0.0) then
        DFACT=SQRT(fk2Rmax)*ACC(N+6)
      else
        DFACT=fk(n)*acc(n+6)
      endif
      IT=DFACT/DFK
      AKEX(N)=0.0
      AFK(N)=DFACT*RMAX
C
C     USE WKB TO FIND APPROX. ABSORPTION
C
      UMXX=DFACT**2
      UMNN=DFACT**2
      T3=0.0
      T4=0.0
      MK=2*kt-1
      I3=MK+2*KMAX*(N-1)
      MK=MK+K-1
      DO 30 M=1,K
      TH=ATAN(U(I3+1)/U(I3  ))*0.5
      IF(TH.GT.0.0) GO TO 28
      IF(U(I3  ).LT.0.0) TH=TH+PI/2.0
   28 CONTINUE
      Y=SQRT(SQRT(U(I3  )**2+U(I3+1)**2))
      X=Y*COS(TH)
      Y=Y*SIN(TH)
      T4=0.5*DRF*(T3+Y)+T4
      FF(MK+2)=T4
      T3=Y
      UMXX= MAX (UMXX,X**2)
      FF(MK+1)=UMXX
      UMNN= MIN (UMNN,U(I3))
      FF(MK  )=UMNN
      I3=I3-2
      MK=MK-3
   30 CONTINUE
      T4=ABS(T4)
      IF(ACC(N+2).NE.0.0) GO TO 32
      T2=(1.0-EXP(-4.0*T4))*ABSC*ACCU
      T3= MIN ( MAX (FUN(UMXX)/DFK,DFACT/DFK+T2)- MIN (FUN(UMNN)/DFK
     1, MAX (0.D0,DFACT/DFK-T2)),DFLOAT(NMAX))-20.0
      ACC(N+2)=1.0
      IF(T3.GT.0.0) ACC(N+2)= MAX (1.0-T3/30.D0,0.D0)
   32 CONTINUE
C
C
      GAM=WIDTH*ACC(N+2)
      M=0
      J=0
      DO 50 I=1,128
      T1=FLOAT(I-IT-1)+1.0e-8
C
C     WAVE NO. SPACING FUNCTION
C
      TEMP=DFK*T1*(1.0-GAM**2/(3.0*(GAM**2+T1**2)))+DFACT
C
      IF(TEMP.gt.DFK*0.5.or.temp.gt.dfact*0.5) then
        J=J+1
        ETEMP(J)=TEMP
        IF(ABS(TEMP-DFACT).LT.0.5*dfk) ITEM(N)=J
        IF(M.EQ.0.AND.TEMP.GT.DFACT) M=J
      endif
   50 CONTINUE
      IF(ITEM(N).EQ.0) ITEM(N)=M
C
C     SKEW WAVE SPREAD FOR MOMENTUM TRANSFER
C
      TEMP=0.0
      T2= MIN (FMA(1),FMA(2))/ MAX (FMA(1),FMA(2))
      IF(FMA(N).LT.FMA(3-N)) TEMP=(FK(N)-T2*FK(3-N))/DFK
      IF(FMA(N).GT.FMA(3-N)) TEMP=(T2*FK(N)-FK(3-N))/DFK
      AUG(1,N)= MAX (ABSC, TEMP+1.0)*ACCU
      AUG(2,N)= MAX (ABSC,-TEMP+1.0)*ACCU
      MMIN=MAX0(ITEM(N)-INT(AUG(1,N)),1)
      MMAX=     ITEM(N)+INT(AUG(2,N))
C
C     FIND MAX AND MIN K**2 AND ABSORPTION AT TURNING POINT
C
      T3=DFACT
      T4=DFACT
      DO 110 LL=1,LPLUS
      L=LPLUS-LL
      I1= MAX (FLOAT(L-3)-ACCL,0.D0)/(DRF*DFACT)
      MK=3*MAX0(MIN0(I1,kt),1)-2
      UMNN=FF(MK  )
      UMXX=FF(MK+1)
      TX=ABS(FF(MK+2))
      T2=1.0-EXP(-4.0*TX)
      T1=DFACT-(T2*DFK)*AUG(1,N)
      T2=DFACT+(T2*DFK)*AUG(2,N)
      IY=ITEM(N)
C
C     SEARCH FOR LOWER LIMIT
C
      T1= MAX (SQRT(FLOAT(L*L+L)/RMAX**2)
     1 , MIN (T1, MIN (FUN(UMNN),T3)))
      I3=IY
      DO 90 I=1,IY
      IF(    T1-ETEMP(I) .GT.0.01) GO TO 90
      I3=I
      GO TO 91
   90 CONTINUE
   91 CONTINUE
      T3= MIN (T3,T1)
      IF(ETA(N).NE.0.0) I3=MIN0(I3,IY-1)
      MMIN=MIN0(MMIN,I3)
      N1(IL+L+1 )=MMIN
      N3(IL+L+1 )=I3
C
C     SEARCH FOR UPPER LIMIT
C
      T2= MAX (T2,FUN(UMXX))
      T2= MAX (T2,T4)
      I4=IY
      DO 95 I=IY,J
      IF(    T2-ETEMP(I) .GT.0.01) GO TO 95
      I4=I
      GO TO 96
   95 CONTINUE
      I4=J
   96 CONTINUE
      T4= MAX (T4,T2)
      IF(ETA(N).NE.0.0) I4=MAX0(I4,IY+1)
      MMAX=MAX0(MMAX,I4)
      N2(IL+L+1 )=MMAX
      N4(IL+L+1 )=I4
C     DIAGNOSTIC PRINTOUT
      IF(ICON(20).EQ.9) WRITE(6,7779)L,I3,I4,T1,T2,TX
 7779 FORMAT(3I6,1P3E12.4)
C
  110 CONTINUE
      TSUM(N)=0.0
      DO 120 M=1,KMAX
      MK=M+M+IZ-1
      U(MK  )=U(MK  )*DR2+1.0
      U(MK+1)=U(MK+1)*DR2
      V(MK  )=V(MK  )*DR2
      V(MK+1)=V(MK+1)*DR2
  120 CONTINUE
      MMIN=max0(MMIN-1,0)
      MMAX=min0(MMAX-MMIN  ,NMAX)
C
C     SET UP MATCHING TABLES
C
      DO 122 M=1,NMAX
      MM=M+IX
      EZER(MM)=0.0
      IF(M.GT.MMAX) GO TO 122
      EZER(MM)=ETEMP(M+MMIN)
      UZ(MM)=DR2*EZER(MM)**2+1.0
      R2(MM)=EZER(MM)*RMAX
      R1(MM)=R2(MM)-EZER(MM)*DRF
      S1(MM)=COS(R1(MM))
      C3(MM)=S1(MM)
      S2(MM)=COS(R2(MM))
      C4(MM)=S2(MM)
      S3(MM)=SIN(R1(MM))
      C1(MM)=-S3(MM)
      S4(MM)=SIN(R2(MM))
      C2(MM)=-S4(MM)
      CALL BES(LPLUS,R1(MM),D1(MM),D3(MM),C1(MM),C3(MM))
      CALL BES(LPLUS,R2(MM),D2(MM),D4(MM),C2(MM),C4(MM))
  122 CONTINUE
C
      ITEM(N)=ITEM(N)-MMIN
      DO 125 LL=1,LPLUS
      N1(IL+LL )=N1(IL+LL )-MMIN
      N2(IL+LL )=N2(IL+LL )-MMIN
      N3(IL+LL )=N3(IL+LL )-MMIN
      N4(IL+LL )=N4(IL+LL )-MMIN
  125 CONTINUE
C
C
      I1=IX+1
      I2=IX+MMAX
      WRITE(6,9200)N,ACC(N),ACC(N+2),ACC(N+4),(EZER(M),M=I1,I2)
      IX=IX+NMAX
      IZ=IZ+KMAX+KMAX
      IL=IL+LMAX
      I1MAX(N)=MMAX
  145 CONTINUE
C
C     L LOOP
C
      I1MIN(1)=1
      I1MIN(2)=1
      DO 1000 LL=1,LPLUS
      L=LL-1
      FL=L
      FLF=(FL**2+FL)/12.0
C
C     SET UP INITIAL CONDITIONS FOR WAVE L
C
      JL=0
      IT=0
      IL=LL
      IE=0
      MAXN=0
C
c     Loop for initial and final wave
c
      DO 900 J=1,2
      KMIN= MAX (FLOAT(LL-IETA3)/(2.0*FK(J)*DRF),1.D0)
      index = mod(kmin-1,5)+1
      MK=2*(KMIN+KMAX*(J-1))-1
      FLM2=FLF/FLOAT(KMIN)**2
      IJ=0
      I3=N3(LL+JL)
      I4=N4(LL+JL)
      JS=NS(J)
      FI=-FS(J)
      DO 150 I=1,JS
      FJ=FL+FI
      ANORM(I)=(WAIT(MK  )**2+WAIT(MK+1)**2)*ETA1**2/WT(INDEX)
      IJ=IJ+1
      FACT(I)=0.5*(FJ**2+FJ-FL**2-FL-FS(J)**2-FS(J))
      IK=I+I-1
      F1(IK  )=0.0
      F1(IK+1)=0.0
      F2(IK  )=ETA1
      IF(FJ.LT.ABS(FL-FS(J))) F2(IK  )=0.0
      F2(IK+1)=0.0
      Q1(IK  )=0.0
      Q1(IK+1)=0.0
      Q2(IK  )=F2(IK  )*(U(MK  )+V(MK  )*FACT(I)-FLM2)
      Q2(IK+1)=F2(IK  )*(U(MK+1)+V(MK+1)*FACT(I)     )
      IF(LL.EQ.2) Q1(IK  )=-F2(IK  )/6.0
      FI=FI+1.0
  150 CONTINUE
      IJ2=IJ+IJ
C
      IZ=1
      T1=WAIT(MK  )*ETA1**2
      T2=WAIT(MK+1)*ETA1**2
      MMIN=N1(LL+JL)
      MMAX=N2(LL+JL)
      MMAX=MIN0(MMAX,NMAX)
      DO 160 N=MMIN,MMAX
      NN=N+MAXN
      B2(NN)=ETA1
      B1(NN)=0.0
      P2(NN)=ETA1*(UZ(NN)-FLM2)
      P1(NN)=0.0
      IF(LL.EQ.2) P1(NN)=-ETA1/6.0
      DO 155 I=1,IJ
      AEX(IZ  )=T1
      AEX(IZ+1)=T2
      IZ=IZ+2
  155 CONTINUE
  160 CONTINUE
C
C     INTEGRATE ALL WAVES FOR GIVEN L
C
      IY=1
      IF(KMIN.GE.kt) GO TO 391
      KMIN=KMIN+1
      DO 380 M=KMIN,kt
      index = mod(m-1,5)+1
      MK=MK+2
      FLM2=FLF/FLOAT(M)**2
      DO 200 I=1,IJ
      IK=I+I-1
      T1=12.*F2(IK  )-10.*Q2(IK  )-Q1(IK  )
      T2=12.*F2(IK+1)-10.*Q2(IK+1)-Q1(IK+1)
      F1(IK  )=F2(IK  )
      F1(IK+1)=F2(IK+1)
      Q3(1)   =    U(MK  )+V(MK  )*FACT(I)-FLM2
      Q3(2)   =    U(MK+1)+V(MK+1)*FACT(I)
      DET=Q3(1)**2+Q3(2)**2
      F2(IK  )=(T1*Q3(1)+T2*Q3(2))/DET
      F2(IK+1)=(T2*Q3(1)-T1*Q3(2))/DET
      FF(IY  )=F2(IK  )*WAIT(MK  )-F2(IK+1)*WAIT(MK+1)
      FF(IY+1)=F2(IK+1)*WAIT(MK  )+F2(IK  )*WAIT(MK+1)
      ANORM(I)=ANORM(I)+(FF(IY  )**2+FF(IY+1)**2)/WT(INDEX)
      Q1(IK  )=Q2(IK  )
      Q1(IK+1)=Q2(IK+1)
      Q2(IK  )=T1
      Q2(IK+1)=T2
      IY=IY+2
  200 CONTINUE
  380 CONTINUE
C
C     ADJUST WEIGHT OF LAST POINT
C
      IY=IY-IJ-IJ
      DO 390 I=1,IJ
      ANORM(I)=ANORM(I)+(FF(IY  )**2+FF(IY+1)**2)/WT(INDEX)
      IY=IY+2
  390 CONTINUE
  391 CONTINUE
C
C     CALCULATE NORM. CONST. AND SCATT. AMPLITUDES
C
      SMAT=0.0
      IZ=LL+LL-1+IE
      ARG=S(IL)-S(IL-LL+1)
      Q(1)=COS(ARG)
      Q(2)=SIN(ARG)
      FJ=FL-FS(J)
      DO 410 I=1,IJ
      IK=I+I-1
      DET=F (IL)*GP(IL)-FP(IL)*G (IL)
      A(1)=(F1(IK  )*GP(IL)-F2(IK  )*G (IL))/DET
      A(2)=(F1(IK+1)*GP(IL)-F2(IK+1)*G (IL))/DET
      B(1)=(F2(IK  )*F (IL)-F1(IK  )*FP(IL))/DET
      B(2)=(F2(IK+1)*F (IL)-F1(IK+1)*FP(IL))/DET
      Q(3)=A(1)+B(2)
      Q(4)=B(1)-A(2)
      DET=Q(3)**2+Q(4)**2
      IF(DET.EQ.0.0) DET=1.0
      T(1)=Q(3)/DET
      T(2)=Q(4)/DET
      CX(IK  )=T(1)*Q(1)-T(2)*Q(2)
      CX(IK+1)=T(1)*Q(2)+T(2)*Q(1)
      ANORM(I)=ANORM(I)*DRF*(CX(IK  )**2+CX(IK+1)**2)
      E(IZ  )=B(1)*T(1)-B(2)*T(2)
      E(IZ+1)=B(1)*T(2)+B(2)*T(1)
      TEMP=E(IZ+1)-E(IZ)**2-E(IZ+1)**2
      SMAT= MAX (SMAT,ABS(4.0*TEMP))
C
C     CALCULATE REACTION SIGMA
C
      if(isym(J) .and. is(J).eq.0) TEMP = 4.0*TEMP 
      TSUM(J)=TSUM(J)+(2.0*FJ+1.0)*TEMP
      IZ=IZ+LPL2
      FJ=FJ+1.0
  410 CONTINUE
C
C     SET EXPANSION RANGE
C
      I6=SMAT*AUG(2,J)+1
      I2=MIN0(NMAX,MAX0(I4,ITEM(J)+I6),I1MAX(J),mmax)
      I1MAX(J)=I2
C
      T1=SQRT(FLOAT(L*L+L)/RMAX**2)
      I1=ITEM(J)-1
      DO 420 N=1,I1
      IF(T1-EZER(N+MAXN).GT.0.01) GO TO 420
      I1=N
      GO TO 421
  420 CONTINUE
  421 CONTINUE
      IF(ETA(J).NE.0.0) I1=MIN0(I1,ITEM(J)-1)
      I5=SMAT*AUG(1,J)
      I1=MAX0(I1,MIN0(ITEM(J)-I5  ,I3),I1MIN(J))
      i1=max0(min0(item(j),i1),1)
      IF(FLOAT(L).GT.AFK(J)) I1=MAX0(ITEM(J),I1MIN(J))
C
C     INTEGRATE BESSEL FUNCTIONS FOR FOURIER-BESSEL EXPANSION
C
      MK=1
      DO 450 M=KMIN,kt
      FLM2=FLF/FLOAT(M)**2
      IZ=1
      DO 440 N=I1,I2
      NN=N+MAXN
      T1=12.*B2(NN)-10.*P2(NN)-P1(NN)
      P1(NN)=P2(NN)
      P2(NN)=T1
      B1(NN)=B2(NN)
      B2(NN)=T1/(UZ(NN)-FLM2)
C
C     COMPUTE SUMS FOR EXPANSION COEFFICIENTS
C
      IY=MK
      DO 425 I=1,IJ
      AEX(IZ  )=AEX(IZ  )+B2(NN)*FF(IY  )
      AEX(IZ+1)=AEX(IZ+1)+B2(NN)*FF(IY+1)
      IZ=IZ+2
      IY=IY+2
  425 CONTINUE
  440 CONTINUE
      MK=MK+IJ2
  450 CONTINUE
C
C     NORMALIZE EXPANSION COEFFICENTS
C
      MMIN=N1(LL+JL)
      MMAX=N2(LL+JL)
      MMAX=MIN0(MMAX,NMAX)
      MMIN=MAX0(I1MIN(J),MMIN)
      MMAX=MIN0(I1MAX(J),MMAX)
      IZ=1
      IAUG=0
      DO 460 N=I1,I2
      NN=N+MAXN
      DET=    S3(NN)*C4(NN)-S4(NN)*C3(NN)
      A1=B1(NN)*C4(NN)-B2(NN)*C3(NN)
      A2=B2(NN)*S3(NN)-B1(NN)*S4(NN)
      FNORM=DRF*DET/(A1*SQRT(1.0+(A2/A1)**2)*EZER(NN))
      IF(ABS(C3(NN)).GT.SCALE) IAUG=N+1
      DO 458 I=1,IJ
      IK=I+I-1
      T1       =(AEX(IZ  )*CX(IK  )-AEX(IZ+1)*CX(IK+1))*FNORM
      T2       =(AEX(IZ  )*CX(IK+1)+AEX(IZ+1)*CX(IK  ))*FNORM
      AEX(IZ  )=T1
      AEX(IZ+1)=T2
      IZ=IZ+2
  458 CONTINUE
  460 CONTINUE
      IF(IAUG.NE.0) then
        I1=IAUG
        I1MIN(J)=I1
      endif
C
C     RECUR BESSEL FUNCTIONS TO NEXT L VALUE
C
C     IRREGULAR FCT AT R=RMAX-DR
C
      FLM=FL+FL+1.0
      DO 500 N=MMIN,MMAX
      NN=N+MAXN
      X1=FLM/R1(NN)
      TX=X1*C3(NN)-C1(NN)
      C1(NN)=C3(NN)
      C3(NN)=TX
C
C     REGULAR FCT AT R=RMAX-DR
C
      A3=D3(NN)
      A2=D1(NN)
      AL=LPLUS+LPLUS+3
      DO 470 M=LL,LPLUS
      AL=AL-2.0
      A1=AL*A2/R1(NN)-A3
      A3=A2
      A2=A1
  470 CONTINUE
      S1(NN)=S3(NN)
      S3(NN)=A3
C
C     IRREGULAR FCT AT R=RMAX
C
      X2=FLM/R2(NN)
      B2(NN)=S2(NN)
      TX=X2*C4(NN)-C2(NN)
      C2(NN)=C4(NN)
      C4(NN)=TX
C
C     REGULAR FCT AT R=RMAX
C
      A3=D4(NN)
      A2=D2(NN)
      AL=LPLUS+LPLUS+3
      DO 475 M=LL,LPLUS
      AL=AL-2.0
      A1=AL*A2/R2(NN)-A3
      A3=A2
      A2=A1
  475 CONTINUE
      S2(NN)=S4(NN)
      S4(NN)=A3
C
  500 CONTINUE
 
C
C     SET UP NORMALIZATION MATRIX AND COEFFICIENT MATRIX
C
      IY=0
      NDIM=I2-I1+1
      INDEX=0
      INDXR=0
      DO 640 N=I1,I2
      NN=N+MAXN
      INDXR=INDEX+1
c	iii1=index+1
      DO 620 M=N,I2
      INDEX=INDEX+1
      IF(M.ne.N) then
        INDXR=INDXR+NDIM
        MM=M+MAXN
C
C       OFF DIAGONAL MATRIX ELEMENT
C
        AMAT(INDEX)=RMAX*(R2(MM)*B2(MM)*S2(NN)-R2(NN)*B2(NN)*S2(MM))
     1      /((R2(NN)**2-R2(MM)**2)*EZER(MM)*EZER(NN))
        AMAT(INDXR)=AMAT(INDEX)
      else
C
C       DIAGONAL MATRIX ELEMENT
C
        AMAT(INDEX)=RMAX*(S2(NN)**2-S4(NN)*B2(NN))/(2.0*EZER(NN)**2)
      endif
  610 CONTINUE
  620 CONTINUE
      INDEX=INDEX+N-I1+1
      IZ=IY
      DO 630 M=1,IJ2
      IZ=IZ+1
      CMAT(IZ)=AEX(IZ)
  630 CONTINUE
      IY=IY+IJ+IJ
  640 CONTINUE
C
C
      DET=DETZ
      IMIN=1
      IERR=0
C
      CALL INVERT(AMAT,CMAT,DET ,NDIM,JS+JS,IMIN,IERR)
C
      AKEX(J)= MAX (AKEX(J),EZER(I2+MAXN))
      I1=I1+IMIN-1
      I1MIN(J)=min0(I1,item(j))
      N1(LL+JL)=I1
      N2(LL+JL)=I2
C
C     CHECK ON CONVERGENCE OF WAVES
C
      IF(ICON(12).ne.0) then
      if(j.eq.1) then
        WRITE(6,'(4h0L =,I4,4X,100(1H*))')L
      endif
      WRITE(6,9201)J,I1,I2,DET,IMIN-1
      IY=(IMIN-1)*IJ2+1
      DO 850 I=1,IJ
      SUM(I)=0.0
      IZ=IY
      DO 840 N=IMIN,NDIM
      TEMP  =       CMAT(IZ     )*AEX(IZ  )+CMAT(IZ+1   )*AEX(IZ+1)
      IF(ANORM(I).NE.0.0) TEST(N)=TEMP/ANORM(I)
      IF(ANORM(I).EQ.0.0) TEST(N)=0.0
      SUM(I)=SUM(I)+TEMP
      IZ=IZ+IJ+IJ
  840 CONTINUE
      IY=IY+2
      IF(ANORM(I).EQ.0.0) then
        SUM(I)=0.0
      else
        SUM(I)=(1.0-SUM(I)/ANORM(I))*100.0
      endif
      IF(ICON(12).ge.2) then
        WRITE(6,9100)ANORM(I),SUM(I),(TEST(N),N=IMIN,NDIM)
      endif
  850 CONTINUE
      IF(ICON(12).EQ.1) WRITE(6,9202)(SUM(I),I=1,IJ)
c      IF(IIP.eq.1.and.icon(12).eq.3) then
      IF(icon(12).eq.3) then
        I4=(IMIN-1)*IJ2
        DO 890 N=IMIN,NDIM
        I3=I4+1
        I4=I4+IJ2
        WRITE(6,9101)I1+N-IMIN,   (CMAT(  I),I=I3,I4)
  890   CONTINUE
      endif
      endif
  899 CONTINUE
C
C     WRITE EXPANSION COEFFICENTS ON TAPE 4
C
      I3=(IMIN-1)*IJ2+1
      I4=I3+(NDIM-IMIN+1)*IJ2-1
      WRITE (4)(CMAT(  I),I=I3,I4)
C
      IL=IL+LPLUS
      IE=IE+NS(J)*LPL2
      MAXN=MAXN+NMAX
      JL=JL+LMAX
  900 CONTINUE
 1000 CONTINUE
C
      IF(ICON(5).ne.0) then
C
C     WRITE ELASTIC PARTIAL WAVE SCATTERING AMPLITUDES
C
C
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      WRITE(6,9600)
      WRITE(6,9601)
      DO 1104 LL=1,LPLUS
      L=LL-1
      I1=L+LL
      I2=I1+LPL2*(NS(1)-1)
      WRITE(6,9602)L, (E(J  ),E(J+1 ),J=I1,I2,LPL2)
      I1=I2+LPL2
      I2=I1+LPL2*(NS(2)-1)
      WRITE(6,9603)L, (E(J  ),E(J+1 ),J=I1,I2,LPL2)
 1104 CONTINUE
      endif
      DO 1150 N=1,2
      TSUM(N)=TSUM(N)*4.0*PI/((2.0*FS(N)+1.0)*FK(N)**2)
 1150 CONTINUE
      WRITE(6,9002)TSUM
      REWIND 4
      if(icon(6).ne.0) then
      IY=1
      IK=1
      DO 2000 J=1,2
      IJ=NS(J)
      T3=S(IK)
      DO 1900 LL=1,LPLUS
      ARG=2.0*(S(IK)-T3)
      Q(3)=COS(ARG)/FK(J)
      Q(4)=SIN(ARG)/FK(J)
      IZ=IY
      DO 1850 I=1,IJ
      T1=Q(3)*E(IZ  )-Q(4)*E(IZ+1)
      T2=Q(3)*E(IZ+1)+Q(4)*E(IZ  )
      if(isym(J) .and. is(J).eq.0) then
        T1 = 2.0*T1
        T2 = 2.0*T2
      endif
      E(IZ  )=T1
      E(IZ+1)=T2
      IZ=IZ+LPL2
 1850 CONTINUE
      IY=IY+2
      IK=IK+1
 1900 CONTINUE
      IY=IZ+2-LPL2
 2000 CONTINUE
      endif
 3000 CONTINUE
c
      RETURN
c
 9002 FORMAT(1H0,9HREACSIG 1,1PE14.4,45X,9HREACSIG 2,1PE14.4)
 9100 FORMAT(6H NORM=,F8.3,6H %ERR=,F8.3,2H/ ,10F9.4/(30X,10F9.4))
 9101 FORMAT(1H0,I3,1P10E12.4/(11X,1P10E12.4))
 9200 FORMAT(26H0WAVE NUMBERS FOR CHANNEL ,I2,2X,9HACC(N  )=,F9.4,9X
     1,9HACC(N+2)=,F9.4,9X,9HACC(N+4)=,F9.4,/(30X,10F9.4))
 9201 FORMAT(6H0CHANN,I2,9H WAVE NO.,I2,1H-,I2,10X,9H     DET=,1PE9.2
     1,9H     AUG=,I2)
 9202 FORMAT(1H+,60X,9H    %ERR=,5F8.3)
 9600 FORMAT(1H0,35HPARTIAL WAVE SCATTERING AMPLITUDES )
 9601 FORMAT(4H   L,20H  REAL D1   IMAG D1 ,20H  REAL D2   IMAG D2 ,20H
     1 REAL D3   IMAG D3 ,4X,4H   L,20H  REAL D1   IMAG D1 ,20H  REAL D2
     2   IMAG D2 ,20H  REAL D3   IMAG D3 )
 9602 FORMAT(1H ,    I3,6F10.6)
 9603 FORMAT(1H+,68X,I3,6F10.6)
 9950 FORMAT(14H0ELAPSED TIME=,F10.3,25H     SUBR. INTEG  ENTERED )
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2),I4,2(1H.,I2))
      END

      SUBROUTINE INVERT(A,C,DET,MAX,IK,L,IERR)
c
      IMPLICIT REAL*8(A-H,O-Z)
C
C     BACKWARDS INVERSION ROUTINE
C
      DIMENSION A(1),C(1)
      TEMP=DET
      DET=1.0
      L=1
C
C     REDUCE DIAGONAL TO 1.0
C
      M1=0
      M2=0
      DO 30 M=1,MAX
      IF(A(M1+M).EQ.0.0) GO TO 21
      B=1.0/A(M1+M)
      M3=M
      DO 15 N=1,MAX
      A(M3  )=A(M3  )*B
      M3=M3+MAX
   15 CONTINUE
      DO 20 K=1,IK
      C(M2+K)=C(M2+K)*B
   20 CONTINUE
   21 CONTINUE
      M1=M1+MAX
      M2=M2+IK
   30 CONTINUE
      IF(MAX.EQ.1) GO TO 201
C
C     REDUCE MATRIX TO LOWER DIAGONAL FORM
C
      M2=MAX*IK+1
      M0=MAX*MAX
      M1=M0
      DO 50 M=2,MAX
      M3=M1-1
      M6=M2
C
      DO 45 K=M,MAX
      B=A(M3  )/A(M1  )
      M4=M3
      M5=M1
      MM=MAX-M+2
      DO 35 N=1,MM
      A(M4  )=A(M4  )-A(M5  )*B
      M4=M4-MAX
      M5=M5-MAX
   35 CONTINUE
      M6=M6-IK
      DO 40 N=1,IK
      C(M6-N)=C(M6-N)-C(M2-N)*B
   40 CONTINUE
      M3=M3-1
   45 CONTINUE
      M1=M1-MAX-1
      M2=M2-IK
   50 CONTINUE
C
C     DETERMINANT CHECK
C
      M1=M0
      K=MAX
      B=DET
      DO 60 M=1,MAX
      DET=DET*A(M1)
      if(det .lt. 0.0) ierr = -1
      IF(DET.LT.TEMP) GO TO 61
      B=DET
      K=K-1
      M1=M1-MAX-1
   60 CONTINUE
   61 CONTINUE
      L=K+1
      DET=B
C
C     SOLVE REMAINING LOWER DIAGONAL
C
      M1=MAX*(L-1)+L
      M2=IK*(L-1)
      DO 100 M=L,MAX
      B=1.0/A(M1)
      DO 70 N=1,IK
      C(M2+N)=C(M2+N)*B
   70 CONTINUE
      M3=M1+1
      M4=M2
      DO 90 N=M,MAX
      IF(N.EQ.M) GO TO 90
      B=A(M3)
      M4=M4+IK
      DO 80 K=1,IK
      C(M4+K)=C(M4+K)-C(M2+K)*B
   80 CONTINUE
      M3=M3+1
   90 CONTINUE
      M1=M1+MAX+1
      M2=M2+IK
  100 CONTINUE
  201 CONTINUE
      RETURN
      END

      SUBROUTINE BES(L,X,F1,F2,G1,G2)
c
      IMPLICIT REAL*8(A-H,O-Z)
      DATA SCALE/1.0E-30/
C
C     COMPUTES REGULAR SPHERICAL BESSEL FUNCTIONS
C     FOR A DEFINITE L AND ARGUMENT
C
      LMAX=2.0*X+10.0
      LMAX=MAX0(LMAX,L+10)
      A3=0.0
      A2=SCALE
      FL=LMAX+LMAX+3
      DO 40 M=1,LMAX
      FL=FL-2.0
      A1=FL*A2/X-A3
      IF(M.NE.LMAX+1-L) GO TO 20
      T1=SCALE/A3
      A3=A3*T1
      A2=A2*T1
      A1=A1*T1
      F2=A3
      F1=A2
      GO TO 30
   20 CONTINUE
      IF(M.GT.LMAX+1-L) GO TO 30
      IF(ABS(A2).LT.1.0) GO TO 30
      A1=A1*SCALE
      A2=A2*SCALE
   30 CONTINUE
      A3=A2
      A2=A1
   40 CONTINUE
      G3=G2/X-G1
      ANORM=1.0/(A2*G3-A3*G2)
      F2=F2*ANORM
      F1=F1*ANORM
      RETURN
      END

