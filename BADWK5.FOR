c$debug
      SUBROUTINE BADWK5
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      COMMON/BDWCK_5/LTR,LTP,JTR,JTP,T(4),QMAX(2),DQ(2)
c
C      CALL KSECOND(TIME)
C      WRITE(6,9950)TIME
C
C
      CALL EXPAND(LTR,DQ,QMAX)
C
C      CALL KSECOND(TIME)
C      WRITE(6,9951)TIME
      CALL XFORM(LTR,LTP,T,DQ)
      RETURN
 9950 FORMAT('0ELAPSED TIME =',F10.3,'     SUBR. EXPAND CALLED')
 9951 FORMAT('0ELAPSED TIME =',F10.3,'     SUBR. XFORM STARTED')
      END

      SUBROUTINE EXPAND(JA,DQ,QMAX)
c
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
      parameter (ioffset = 1025, iorg_ff = 5)
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical interp
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
c
      dimension ub(ispc0), ff(2*ioffset)
     1,         qmax(2),wt(5),s1(2),s2(2),c1(2),c2(2),ja(2),dq(2)
     2,         delk(2), Df(10) 
c
      equivalence (Space0(   1),ub(1)), (Space2(   1),ff(1))
c
      DATA ETA1, pi/1.0E-08, 3.14159265359E0/
C
      SQR4PI=SQRT(4.*PI)
      WT(1)=375./288.
      WT(4)=WT(1)
      WT(2)=250./288.
      WT(3)=WT(2)
      WT(5)=190./288.
      IORG = iorg_ff
C
C     FORM FACTOR LOOP
C
      DO 200 N=1,2
      DRX=DR(N+2)
      DK=0.05
   25 DELK(N)=DK
      IMAX=QMAX(N)/DK+5.0
      IF(IMAX.LE.505) GO TO 30
      DK=QMAX(N)/500.
      WRITE(6,9104)DK
      GO TO 25
   30 IF(N.EQ.1) GO TO 36
      IF(ICON(11).NE.0) GO TO 101
   36 CONTINUE
      READ (9)KT,FN,ETA(4),FK(4),(UB(I),I=1,KT), Df 
      KT=5*(KT/10)
      LTR=JA(N)
      FLF=LTR*(LTR+1)
C
      R=0.0
      DO 37 M=1,KT
      R=R+DRX
      UB(M+M  )=UB(M+M-1)*R
   37 CONTINUE
      interp=.true.
c	write(20,'(a,2i4,1p4e12.4)')' N, KT, UB ',N,KT
c     1        ,(UB(I),I=2*kt-4,2*kt,2) 
      FKX2 = ((UB(2*KT-4)+UB(2*KT  )-2.0*UB(2*KT-2))/UB(2*KT-2)
     1   -FLF/FLOAT(KT-1)**2)/DRX**2
      FKX = SQRT(ABS(FKX2))
c	write(6,'(a,i4,1p4e14.6)')' N, FKX, FKX2 ',N,FKX,FKX2
C
      R1=FKX*FLOAT(KT)*DRX
      PHAS=PHASEF(LTR)
      UBZ=0.0
      H1=1.0
      H2=1.0
      IF(LTR.ne.0) then
        DO 38 L=1,LTR
        FACT=L+L-1
        T3=FACT*H2/R1+H1
        H1=H2
        H2=T3
   38   CONTINUE
      endif
      H3=FLOAT(LTR)/R1+H1/H2
C
      Q=0.0
      IF(LTR.EQ.1) then
        ETA2=-ETA1/6.0
      else
        ETA2=0.0
      endif
      FLFX=FLF/12.
C
C     MOMENTUM TRANSFER LOOP
C
      DO 100 I=1,IMAX
      q=q+dk
      IF(Q*DRX.gt.0.25.and.kt.le.1600) then
c       Spread out table
        DO 41 M=1,KT
        MM=KT-M
        UB(MM+MM+2)=UB(MM+1)
   41   CONTINUE
        interp=.true.
      endif
      if(interp) then
        interp=.false.
        DRX=DRX*0.5
        UB(KT+KT+2)=UB(KT+KT)**2/UB(KT+KT-2)
c       Interpolate at in between points
        DO 43 M=3,KT
        UB(M+M-1)=(9.*(UB(M+M-2)+UB(M+M  ))-UB(M+M-4)-UB(M+M+2))*.0625
   43   CONTINUE
        UB(1)=(9.*(UBZ+UB(2))-UB(2)*PHAS-UB(4))*.0625
        UB(3)=(9.*(UB(2)+UB(4))-UBZ-UB(6))*.0625
        KT=2*KT
      endif
      SUM=0.0
      UZ=1.0+(Q*DRX)**2/12.
C
C     SET INITIAL CONDITIONS FOR BESSEL FCTS
C
      B1=0.0
      B2=ETA1
      P1=ETA2
      P2=ETA1*(UZ-FLFX)
      R=0.0
      INDEX=1
      SUM=B2*UB(1)*WT(1)
C
C     INTEGRATION LOOP
C
      DO 50 M=2,KT
      T1=12.*B2-10.*P2-P1
      P1=P2
      P2=T1
      B1=B2
      B2=T1/(UZ-FLFX/FLOAT(M)**2)
      INDEX=INDEX+1
      SUM=SUM+B2*UB(M)*WT(INDEX)
      IF(INDEX.EQ.5) INDEX=0
   50 CONTINUE
      SUM=SUM-B2*UB(KT)*WT(5)*0.5
C
C     BESSEL FCTS FOR L.EQ.-1 AND L.EQ.0
C
      R2=Q*DRX*FLOAT(KT)
      R1=R2-Q*DRX
      S1(1)=COS(R1)
      S1(2)=SIN(R1)
      C1(1)=-S1(2)
      C1(2)= S1(1)
      S2(1)=COS(R2)
      S2(2)=SIN(R2)
      C2(1)=-S2(2)
      C2(2)= S2(1)
      IF(LTR.ne.0) then
C
C     RECUR BESSEL FCTS FOR L.GT.0
C
        DO 60 L=1,LTR
        FACT=L+L-1
        T1=FACT*S1(2)/R1-S1(1)
        S1(1)=S1(2)
        S1(2)=T1
        T1=FACT*C1(2)/R1-C1(1)
        C1(1)=C1(2)
        C1(2)=T1
        T2=FACT*S2(2)/R2-S2(1)
        S2(1)=S2(2)
        S2(2)=T2
        T2=FACT*C2(2)/R2-C2(1)
        C2(1)=C2(2)
        C2(2)=T2
   60   CONTINUE
      endif
      S3=S2(1)-FLOAT(LTR)*S2(2)/R2
C
C     NORMALIZE BESSEL FCTS
C
      DET = S2(2)*C1(2)-S1(2)*C2(2)
      A1 =  B1*C2(2)-B2*C1(2)
      A2 =  B2*S1(2)-B1*S2(2)
      FNORM=DRX*DET/SQRT(A1**2+A2**2)
      SUM=SUM*FNORM
C
C     ADD ON FROM RMAX TO INFINITY
C
      SUM=SUM+UB(KT)*(Q*S3+FKX*H3*S2(2))/(Q**2+FKX2)
      IN=IORG+I+I
c     Point at dq stored at in = 5
      FF(IN)=SUM*SQR4PI/Q**(LTR+1)
  100 CONTINUE
      GO TO 106
C
C     READ IN EXTERNAL FOURIER TRANSFORM
C
  101 CONTINUE
      READ (9)KT,FN,ETA(4),FK(4),(UB(I),I=1,KT)
      DK=FN
      DELK(N)=DK
      IMAX=KT
      IMAX=MIN0(IMAX,510)
      DO 105 M=1,IMAX
      MK=IORG+2*M
c     Point at dq stored at mk = 5
      FF(MK)=UB(M)
  105 CONTINUE
  106 CONTINUE
C
c     Point at    dq stored at in = 7
c     Point at     0 stored at in = 5
c     Point at -  dq stored at in = 3
c     Point at -2*dq stored at in = 1
      FF(IORG  )=(4.*FF(IORG+2)-FF(IORG+4))/3.
      FF(IORG-2)=FF(IORG+2)
      FF(IORG-4)=FF(IORG+4)
      WRITE(6,9102)N, LTR 
      I1=IORG
      Q=0.0
      DO 110 I=1,IMAX,10
      I2=MIN0(I1+19,IORG+2*IMAX)
      WRITE(6,9101)Q,(FF(M),M=I1,I2,2)
      I1=I2+1
      Q=Q+10.*DK
  110 CONTINUE
      IMAX=IMAX-3
C
C     6 POINT LAGRANGIAN INTERPOLATION LOOP
C
      DO 120 I=1,IMAX
      IN=IORG+I+I-1
      FF(IN)=(450.*(FF(IN+1)+FF(IN-1))-75.*(FF(IN+3)+FF(IN-3))
     1+9.*(FF(IN+5)+FF(IN-5)))/768.
  120 CONTINUE
      FF(IORG-1)=FF(IORG+1)
      IORG=IORG+ioffset
  200 CONTINUE
      DQ(1)=DELK(1)/2.0
      DQ(2)=DELK(2)/2.0
      RETURN
c
 9101 FORMAT(1H ,F8.2,1P10E12.4)
 9102 FORMAT('0    Q    FOURIER TRANSFORM FOR FORM FACTOR ',I3
     1    ,'    L = ', i3)
 9104 FORMAT('0TOO LARGE MOM. RANGE FOR FOURIER XFORM, INCREASE DK = '
     1,F8.4,8X,50(1H$))
      END

      SUBROUTINE XFORM(LTR,LTP,TX,DK)
c
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
c
      IMPLICIT REAL*8(A-H,O-Z)
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
      Common/n1n2ez/n1(200,2),n2(200,2),ezer(100)
c
      dimension stor1(imat*imat,2), stor2(imat*imat,2)
     1,         f(1025),g(1025),xg(250),wg(250),p(250,4),e(imat,8)
      DIMENSION TX(4),T(4),DQ(2),DK(2)
c
      equivalence (Space0(   1),stor1(1,1)), (Space1(   1),stor2(1,1))
     1,           (Space2(   1),f(1)),       (Space2(1026),g(1))
     2,           (Space3(   1),p(1,1)),     (Space3(1001),xg(1))
     3,           (Space3(1251),wg(1)),      (Space3(1501),e(1,1))
c
      DATA PI/3.14159265359/
C
C
      LPLUSX=LPLUS+LTR+LTP
      DQ(1)=DK(1)
      DQ(2)=DK(2)
      F2=0.0
      G2=0.0
      FACTOR=1.0/SQRT(16.*PI)
      DO 10 I=1,4
      T(I)=TX(I)
   10 CONTINUE
      M1=0
      DO 80 J=1,2
      DO 70 N=1,NMAX
      NN=N+M1
      IF(EZER(NN).EQ.0.0) GO TO 75
      E(N,J+4)=EZER(NN)*T(J  )
      E(N,J  )=E(N,J+4)**2
      E(N,J+6)=EZER(NN)*T(J+2)
      E(N,J+2)=E(N,J+6)**2
   70 CONTINUE
   75 CONTINUE
      M1=M1+NMAX
   80 CONTINUE
      NGAUS=MIN0((LPLUSX+10*(ICON(19)+1))/2,250)
      LGAUS=NGAUS+NGAUS
      CALL LEGAUS(LGAUS,XG,WG)
      DO 100 I=1,NGAUS
      P(I,1)=WG(I)
      P(I,2)=WG(I)*XG(I)
      P(I,3)=(3.0*XG(I)*P(I,2)-    P(I,1))/2.0
      P(I,4)=(5.0*XG(I)*P(I,3)-2.0*P(I,2))/3.0
  100 CONTINUE
      DO 500 LL=1,LPLUSX,4
      IF(LL.GT.LPLUS) GO TO 121
      NIA=N1(LL,1)
      NIB=N2(LL,1)
      NFA=N1(LL,2)
      NFB=N2(LL,2)
      LTX=LTR+LTP+2
      DO 120 I=1,LTX
      INDA=LL-I+2
      IF(INDA.LT.1) GO TO 115
      IF(INDA.GT.LPLUS) GO TO 115
      NIA=MIN0(NIA,N1(INDA,1))
      NFA=MIN0(NFA,N1(INDA,2))
      NIB=MAX0(NIB,N2(INDA,1))
      NFB=MAX0(NFB,N2(INDA,2))
  115 CONTINUE
      INDA=LL+I+1
      IF(INDA.GT.LPLUS) GO TO 120
      NIA=MIN0(NIA,N1(INDA,1))
      NFA=MIN0(NFA,N1(INDA,2))
      NIB=MAX0(NIB,N2(INDA,1))
      NFB=MAX0(NFB,N2(INDA,2))
  120 CONTINUE
  121 CONTINUE
      INDEX=0
      DO 400 NI=NIA,NIB
      DO 350 NF=NFA,NFB
      T1=E(NI,1)+E(NF,2)
      T3=E(NI,3)+E(NF,4)
      T5=2.0*E(NI,5)*E(NF,6)
      T6=2.0*E(NI,7)*E(NF,8)
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      SUM4=0.0
      DO 250 I=1,NGAUS
      T2=T5*XG(I)
      T4=T6*XG(I)
C
C     4 POINT LAGRANGE INTERPOLATION
C
c     Angles .le. 90 degrees
c
      D1=SQRT(T1-T2)/DQ(1)
      M1=D1
      D1=D1-FLOAT(M1)
      F1=-D1*(D1-1.)*(D1-2.)*F(M1+4)/6.+(D1**2-1.)*(D1-2.)*F(M1+5)/2.
     1-D1*(D1+1.)*(D1-2.)*F(M1+6)/2.+D1*(D1**2-1.)*F(M1+7)/6.
C
      D1=SQRT(T3-T4)/DQ(2)
      M1=D1
      D1=D1-FLOAT(M1)
      G1=-D1*(D1-1.)*(D1-2.)*G(M1+4)/6.+(D1**2-1.)*(D1-2.)*G(M1+5)/2.
     1-D1*(D1+1.)*(D1-2.)*G(M1+6)/2.+D1*(D1**2-1.)*G(M1+7)/6.
C
c     Angles .ge. 90 degrees
c
      IF(ICON(13).NE.0) GO TO 160
      D1=SQRT(T1+T2)/DQ(1)
      M1=D1
      D1=D1-FLOAT(M1)
      F2=-D1*(D1-1.)*(D1-2.)*F(M1+4)/6.+(D1**2-1.)*(D1-2.)*F(M1+5)/2.
     1-D1*(D1+1.)*(D1-2.)*F(M1+6)/2.+D1*(D1**2-1.)*F(M1+7)/6.
C
      D1=SQRT(T3+T4)/DQ(2)
      M1=D1
      D1=D1-FLOAT(M1)
      G2=-D1*(D1-1.)*(D1-2.)*G(M1+4)/6.+(D1**2-1.)*(D1-2.)*G(M1+5)/2.
     1-D1*(D1+1.)*(D1-2.)*G(M1+6)/2.+D1*(D1**2-1.)*G(M1+7)/6.
C
  160 CONTINUE
      SUM1=SUM1+(F1*G1+F2*G2)*P(I,1)
      SUM3=SUM3+(F1*G1+F2*G2)*P(I,3)
      SUM2=SUM2+(F1*G1-F2*G2)*P(I,2)
      SUM4=SUM4+(F1*G1-F2*G2)*P(I,4)
  250 CONTINUE
      INDEX=INDEX+1
      STOR1(INDEX,1)=SUM1*FACTOR
      STOR1(INDEX,2)=SUM2*FACTOR
      STOR2(INDEX,1)=SUM3*FACTOR
      STOR2(INDEX,2)=SUM4*FACTOR
  350 CONTINUE
  400 CONTINUE
      FL1=LL+2
      DO 421 J=1,4
      FL=FL1+1.0
      FL2=FL1+FL
      DO 420 I=1,NGAUS
      P(I,1)=P(I,2)
      P(I,2)=P(I,3)
      P(I,3)=P(I,4)
      P(I,4)=(FL2*XG(I)*P(I,3)-FL1*P(I,2))/FL
  420 CONTINUE
      FL1=FL1+1.0
  421 CONTINUE
      WRITE(3)INDEX,NIA,NIB,NFA,NFB,(STOR1(I,1),I=1,INDEX)
      IF(LL+1.GT.LPLUSX) GO TO 500
      WRITE(3)INDEX,NIA,NIB,NFA,NFB,(STOR1(I,2),I=1,INDEX)
      IF(LL+2.GT.LPLUSX) GO TO 500
      WRITE(3)INDEX,NIA,NIB,NFA,NFB,(STOR2(I,1),I=1,INDEX)
      IF(LL+3.GT.LPLUSX) GO TO 500
      WRITE(3)INDEX,NIA,NIB,NFA,NFB,(STOR2(I,2),I=1,INDEX)
  500 CONTINUE
      REWIND 3
      RETURN
      END
