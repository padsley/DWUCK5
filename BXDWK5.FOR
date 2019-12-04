
      SUBROUTINE BXDWK5
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
c
      IMPLICIT REAL*8(A-H,O-Z)
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
      Common/bdwck_5/ ltr,ltp,jtr,jtp,t(4),qmax(2),dp,pn,dq,qn
c
      dimension si(6400),ub(1600),f1(1025),f2(1025),xq(100),wq(100)
     1,         ff(6400),iqn(5),qnum(4,2),q(4),con(4),zz(2),cf(5)
c
      equivalence (Space1(   1),si(1)), (Space0(   1),ub(1))
     1,           (Space0(1601),f1(1)), (Space0(2626),f2(1))
     2,           (Space0(4001),xq(1)), (Space0(4401),wq(1))
     3,           (Space2(   1),ff(1))
c
c     SI   is target          form factor with dim KSI
c     FF   is project.*target form factor with dim KFF
c     F1   is single particle form factor with dim KF2
c
      data ksi, kff, kf2/ispc0, ispc2, 1020/
c
C      CALL KSECOND(TIME)
C      WRITE(6,9950)TIME, 'Bxdwk5'
      CF(1)=375./288.
      CF(4)=CF(1)
      CF(2)=250./288.
      CF(3)=CF(2)
      CF(5)=190./288.
      SCALE=ABS(FK(1)-FK(2))
C
C     TEST FOR OPTIONS FOR READ IN OF STORED FORM FACTORS
      IF(ICON(11).EQ.2)  GO TO 500
      IF(ICON(11).EQ.3) GO TO 400
      IQN(1)=LTR
      IQN(2)=IS(3)
      IQN(3)=JTR
      IQN(5)=0
C     X IS SINGLE PARTICLE MOMENTUM COORDINATE
      QX=07.5
      DX=0.0375
      KX= MAX (QMAX(1)/2.0,QX)/DX+1.E-6
      KX=MIN0(KX,KF2/2)
      DX=QX/FLOAT(KX)
C
C     DQ AND QN  ARE PROJECTILE RELATIVE COORDINATE SPACING AND NUMBER
C     DP AND PN  ARE PROJ CENTER OF MASS COORDINATE SPACING AND NUMBER
C     RELATIVE MOMENTUM POINTS AND WEIGHTS
      KQ=QN
      I=0
      DO 250 MQ=1,KQ
      I=I+1
      XQ(MQ)=DQ*FLOAT(MQ)
      WQ(MQ)=CF(I)*DQ*XQ(MQ)**2
      IF(I.EQ.5) I=0
  250 CONTINUE
C     COMPUTE TARGET CENTER OF MASS INTERVAL AND NUMBER OF POINTS
      DQRMIN=.02
      KPT=PN
      KP=QMAX(2)/DP+4.0
      KP=MIN0(KP,KPT)
      KR=MIN0(KSI/KQ,KFF/(KP+1)-1)
      DQR= MAX (QMAX(1)/FLOAT(KR-2),DQRMIN)
      KR=QMAX(1)/DQR+2.0001
      PN=KP
      RN=KR
      WRITE(6,9101)
      WRITE(6,9100)DP,PN,DQ,QN,DQR,RN,DX,FLOAT(KX)
C
      DO 40 M=1,KSI
      SI(M)=0.0
   40 CONTINUE
   50 CONTINUE
C
      DO 200 IC=1,2
      IF(IC.EQ.1) WRITE(6,9998)
C     Read in single particle wave functions in r-space 
      READ (9) Q,ZZ(2),TSPIN,CNTRL,OPT,DRF,KT,(UB(2*I-1),I=1,KT)
      IF(CNTRL.EQ.0.0) GO TO 301
      DO 80 I=1,4
      QNUM(I,2)=Q(I)
   80 CONTINUE
      ICN=ABS(CNTRL)
      L=Q(2)
c     Compute fourier transform of single particle wave functions
      CALL EXPAN2(KT,KMAX,L,DRF,DX,UB,F2,KX,ICON(4))
      IF(IC.EQ.2) GO TO 200
      KX2=KX+KX
      DO 100 M=1,KX2
      F1(M-1)=F2(M-1)
  100 CONTINUE
      DO 110 M=1,4
      QNUM(M,1)=QNUM(M,2)
  110 CONTINUE
      ZZ(1)=ZZ(2)
      IF(ICN.EQ.1) GO TO 201
  200 CONTINUE
  201 CONTINUE
      IF(TSPIN.NE.0.0) GO TO 205
      IF(ZZ(1).EQ.ZZ(2)) TSPIN=1.0
  205 IQN(4)=TSPIN
      CON(1)=DQR
      CON(2)=KR
      CON(3)=QN
C
C     TEMP=PHASEF(LTP+IS(3)/2)
C    1    /SQRT(FLOAT((2*LTP+1)*(  IS(3)+1)*(2*LTR+1)))
C    1    /RACAH(2*LTP,2*LTR,JTP,JTR,2*LTR,IS(3))
C     OPT=OPT*TEMP
C     WRITE(6,9952)TEMP
C9952 FORMAT(15H0SCALE FACTOR = ,1PE13.4)
c
c     Compute tnt target form factor
c
      SCALE=ABS(FK(1)-FK(2))
      CALL DSTRI2(IQN,0.5*DX,2*KX,F1,SI,QNUM,OPT,XQ,CON,SX,SCALE)
      IF(CNTRL.GT.0.0) GO TO 50
  301 CONTINUE
      GO TO 440
C
C     READ IN TARGET FORM FACTOR STORED ON DISK
C
  400 CONTINUE
C     DQ  AND QN:  TARGET CENTER OF MASS
C     DQR AND RN:  TARGET RELATIVE MOTION
      READ (8)DQ,DQR,QN,RN
      KR=RN
      KQ=QN
      KF=KR*KQ
      READ (8)(SI(M),M=1,KF)
      WT=2.0
      DO 420 MQ=1,KQ
      WT=6.0-WT
      X=DQ*FLOAT(MQ)
      WQ(MQ)=WT*(DQ/3.0)*X**2
  420 CONTINUE
  440 CONTINUE
C
C     WRITE OUT TARGET FORM FACTOR
C
      IF(ICON(18).GE.2) THEN
      WRITE(6,9041)
      KC=ICON(18)
      QR=DQR
      IA=1
      DO 445 N=1,KR,KC
      IB=IA+KQ-1
      WRITE(6,9051)QR,(SI(M),M=IA,IB)
      QR=QR+DQR*FLOAT(KC)
      IA=IB+(KC-1)*KQ+1
  445 CONTINUE
      ENDIF
C
      RN=KR
C
C     WRITE TARGET FORM FACTOR ON DISK
C
      IF(ICON(20).EQ.8) THEN
        WRITE(8)DQ,CON(1),QN,CON(2)
        KF=KR*KQ
        WRITE(8)(SI(M),M=1,KF)
      ENDIF
C
      READ (9)DQ,DP,QN,PN
C
C     MULTIPLY PROJECTILE AND TARGET FORM FACTORS
C     FF(P,Q) = SUM OVER q [PROJ(P,q) * TARG(q,Q)]
C     RANGES: P IS KP  PROJECTILE CENTER OF MASS
C             q IS KQ  RELATIVE MOTION
C             Q IS KR  TARGET     CENTER OF MASS
C
      J=KR+1
      DO 460 N=1,KPT
      J=J+1
C     MAKE SURE ALL INPUT IS READ
C     READ IN PROJECTILE FORM FACTOR ONE ROW AT A TIME
      READ (9)(XQ(M),M=1,KQ)
      IF(N.LE.KP) THEN
        IB=0
        IA=J
c	  an=si(1)/exp(-2.89*(dqr**2/4.0+dq**2))
        DO 450 M=1,KR
        J=J+1
        FF(J)=0.0
c	  st=an*exp(-2.89*(dqr*float(m))**2/4.0)
        DO 449 MR=1,KQ
        IB=IB+1
c	  si(ib)=st*exp(-2.89*(dq*float(mr))**2)
        FF(J)=FF(J)+SI(IB)*WQ(MR)*XQ(MR)
  449   CONTINUE
  450   CONTINUE
        FF(IA)=(4.0*FF(IA+1)-FF(IA+2))/3.0
      ENDIF
  460 CONTINUE
C
      KR=KR+1
      KP=KP+1
      DO 480 M=1,KR
      FF(M)=(4.0*FF(M+KR)-FF(M+KR+KR))/3.0
  480 CONTINUE
C
C     SAVE TARGET FORM FACTOR ON DISK
C
      IF(ICON(11)/2.NE.2) GO TO 600
      WRITE(8)KR,KP,DQR,DP
      KRKP=KR*KP
      WRITE(8)(FF(M),M=1,KRKP)
      GO TO 600
  500 CONTINUE
      READ (8)KR,KP,DQR,DP
      KRKP=KR*KP
      READ (8)(FF(M),M=1,KRKP)
  600 CONTINUE
C
C     WRITE OUT FORM FACTOR PRODUCT
C
      IF(ICON(18).EQ.1.OR.ICON(18).EQ.3) THEN
      WRITE(6,9040)
      P=0.0
      IA=1
      DO 690 N=1,KP
      IB=IA+KR-1
      WRITE(6,9050)P,(FF(M),M=IA,IB)
      IA=IB+1
      P=P+DP
  690 CONTINUE
      ENDIF
C
C      CALL KSECOND(TIME)
C      WRITE(6,9950)TIME, 'Xform2'
C
      CALL XFORM2(LTR,LTP,T,DQR,DP,KR,KP)
C
      RETURN
 9040 FORMAT('1   Q ====>','  FORM FACTOR PRODUCT')
 9041 FORMAT('1   q ====>','  TARGET FORM FACTOR')
 9050 FORMAT(4H P =,F6.2,1P10E12.4/(10X,1P10E12.4))
 9051 FORMAT(4H QR=,F6.2,1P10E12.4/(10X,1P10E12.4))
 9100 FORMAT(8H    DP =,F7.5,3X,4HNP =,F6.0,4X,4HDq =,F7.5
     1,       3X,4HNq =,F6.0,4X,4HDQ =,F7.5,3X,4HNQ =,F6.0
     2,       4X,4HDx =,F7.5,3X,4HKx =,F6.0)
 9101 FORMAT('0   Projectile c.m. momenta '
     1,'    Proj.-target rel. momenta','   Target c.m. momenta     '
     2,'    Single particle   momenta ')
 9950 FORMAT(14H0ELAPSED TIME=,F10.3,11H     SUBR. ,A5,8H STARTED )
 9998 FORMAT(1H )
      END

      SUBROUTINE EXPAN2(KT,KMAX,LTR,DR,DK,UB,FF,KFX,ICON4)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      DIMENSION WT(5),S1(2),S2(2),C1(2),C2(2),FF(2048)
     1         ,UB(400),C(2,2),S(2,2),RX(2)
c
      EQUIVALENCE (S1(1),S(1,1)),(S2(1),S(1,2)),(R1,RX(1)),(R2,RX(2))
     1  ,(C1(1),C(1,1)),(C2(1),C(1,2))
c
      DATA ETA1, SQR4PI/1.0e-15, 3.54490770D0/
      DATA NCNTR/0/
C
      WT(1)=375./288.
      WT(2)=250./288.
      WT(3)=250./288.
      WT(4)=375./288.
      WT(5)=190./288.
      NCNTR=NCNTR+1
      IORG=0
      DRX=DR
C
C     FORM FACTOR LOOP
C
      IMAX=KFX+3
      KT=5*(KT/10)
      FLF=LTR*(LTR+1)
C
      R=0.0
      DO 37 M=1,KT
      R=R+DRX
      UB(M+M  )=UB(M+M-1)*R
   37 CONTINUE
      FKX=SQRT(ABS((UB(KT+KT-4)+UB(KT+KT  )-2.0*UB(KT+KT-2))/UB(KT+KT-2)
     1-FLF/FLOAT(KT-1)**2))/DRX
C
      R1=FKX*FLOAT(KT)*DRX
      PHAS=PHASEF(LTR)
      UBZ=0.0
      IF(LTR.EQ.0) UBZ=(4.*UB(2)-UB(4))/3.
      H1=1.0
      H2=1.0
      IF(LTR.EQ.0) GO TO 39
      DO 38 L=1,LTR
      FACT=L+L-1
      T3=FACT*H2/R1+H1
      H1=H2
      H2=T3
   38 CONTINUE
   39 CONTINUE
      H3=FLOAT(LTR)/R1+H1/H2
C
      Q=0.0
      ETA2=0.0
      IF(LTR.EQ.1) ETA2=-ETA1/6.0
      FLFX=FLF/12.
C
C     MOMENTUM TRANSFER LOOP
C
 
      DO 100 I=1,IMAX
      Q=Q+DK
      IF(I.EQ.1) GO TO 42
      IF(Q*DRX.LT.0.25) GO TO 44
      IF(KT.GT.KMAX) GO TO 44
      DO 41 M=1,KT
      MM=KT-M
      UB(MM+MM+2)=UB(MM+1)
   41 CONTINUE
   42 CONTINUE
      DRX=DRX*0.5
      UB(KT+KT+2)=UB(KT+KT)**2/UB(KT+KT-2)
      DO 43 M=3,KT
      UB(M+M-1)=(9.*(UB(M+M-2)+UB(M+M  ))-UB(M+M-4)-UB(M+M+2))*.0625
   43 CONTINUE
      UB(1)=(9.*(UBZ+UB(2))-UB(2)*PHAS-UB(4))*.0625
      UB(3)=(9.*(UB(2)+UB(4))-UBZ-UB(6))*.0625
      KT=KT+KT
   44 CONTINUE
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
      IF(LTR.EQ.0) GO TO 85
C
C     RECUR BESSEL FCTS FOR L.GT.0
C
      IF(R1.GT.FLOAT(LTR)) THEN
      DO 60 L=1,LTR
      FACT=L+L-1
      T1=FACT*S1(2)/R1-S1(1)
      S1(1)=S1(2)
      S1(2)=T1
      T2=FACT*S2(2)/R2-S2(1)
      S2(1)=S2(2)
      S2(2)=T2
   60 CONTINUE
C
      ELSE
      LMAX=2.0*R1+10.
      LMAX=MAX0(LMAX,LTR+10)
      DO 75 IS=1,2
      X=RX(IS)
      A3=0.0
      A2=1.0E-30
      FL=LMAX+LMAX+3
      DO 70 M=1,LMAX
      FL=FL-2.0
      A1=FL*A2/X-A3
      IF(M.EQ.LMAX+1-LTR) THEN
      F2=A2
      F1=A1
      ENDIF
      A3=A2
      A2=A1
   70 CONTINUE
      T3=C(2,IS)/X-C(1,IS)
      T2=C(2,IS)
      T1=1.0/(A2*T3-A3*T2)
      S(2,IS)=F2*T1
      S(1,IS)=F1*T1
   75 CONTINUE
      ENDIF
C
      DO 80 L=1,LTR
      FACT=L+L-1
      T1=FACT*C1(2)/R1-C1(1)
      C1(1)=C1(2)
      C1(2)=T1
      T2=FACT*C2(2)/R2-C2(1)
      C2(1)=C2(2)
      C2(2)=T2
   80 CONTINUE
C
   85 CONTINUE
      S3=S2(1)-FLOAT(LTR)*S2(2)/R2
C
C     NORMALIZE BESSEL FCTS
C
      DET=ABS(S1(2)*C2(2)-S2(2)*C1(2))
      A1 =B1*C2(2)-B2*C1(2)
      A2 =B2*S1(2)-B1*S2(2)
      FNORM=DRX*DET/SQRT(A1**2+A2**2)
      SUM=SUM*FNORM
C
C     ADD ON FROM RMAX TO INFINITY
C
      SUM=SUM+UB(KT)*(Q*S3+FKX*H3*S2(2))/(Q**2+FKX**2)
      IN=IORG+I+I
      FF(IN)=SUM*SQR4PI/Q
C	ff(in)=(1.7*sqr4pi)**1.5*exp(-(1.7*q)**2/2.0)
   95 CONTINUE
  100 CONTINUE
C
      IF(LTR.EQ.0) THEN
        FF(IORG  )=(4.*FF(IORG+2)-FF(IORG+4))/3.
      ELSE
        FF(IORG)=0.0
      ENDIF
      WRITE(6,9102)1,NCNTR
      IF(ICON4.EQ.2) WRITE(6,9103) FF(IORG+2)/DK**LTR
      IF(ICON4.EQ.2) GO TO 111
      I1=0
      Q=0.0
      DO 110 I=1,IMAX,10
      I2=MIN0(I1+19,IORG+IMAX+IMAX)
      WRITE(6,9101)Q,(FF(M),M=I1,I2,2)
      I1=I2+1
      Q=Q+10.*DK
  110 CONTINUE
  111 CONTINUE
      IMAX=IMAX-3
      FF(IORG-2)=FF(IORG+2)*PHASEF(LTR)
      FF(IORG-4)=FF(IORG+4)*PHASEF(LTR)
C
C     6 POINT LAGRANGIAN INTERPOLATION LOOP
C
      DO 120 I=1,IMAX
      IN=IORG+I+I-1
      FF(IN)=(450.*(FF(IN+1)+FF(IN-1))-75.*(FF(IN+3)+FF(IN-3))
     1+9.*(FF(IN+5)+FF(IN-5)))/768.
  120 CONTINUE
      RETURN
 9101 FORMAT(1H ,F8.2,1P10E11.3)
 9102 FORMAT(43H     Q    FOURIER TRANSFORM FOR FORM FACTOR,I3,1H.,I3)
 9103 FORMAT('+',60X,'FT(0)=',1PE10.3)
      END

      SUBROUTINE XFORM2(LTR,LTP,T,DK1,DK2,KRR,KP)
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
c
      IMPLICIT REAL*8(A-H,O-Z)
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
      Common/n1n2ez/n1(nls,2),n2(nls,2),ezer(100)
c
      dimension stor1(imat*imat,2) ,stor2(imat*imat,2)
     1,         ff(ispc2),xg(250),wg(250),p(250,4),e(imat,8)
     2,         t(4),dq(2),fx(4),fy(4)
c
      equivalence (Space0(   1),stor1(1,1)), (Space1(   1),stor2(1,1))
     2,           (Space2(   1),ff(1))
     3,           (Space3(   1),p(1,1)),     (Space3(1001),xg(1))
     4,           (Space3(1251),wg(1)),      (Space3(1501),e(1,1))
c
      DATA NPLM, pi/250, 3.14159265359/
C
      LPLUSX=LPLUS+LTR+LTP
      DQ(1)=DK1
      DQ(2)=DK2
      F2=0.0
      G2=0.0
      FACTOR=1.0/SQRT(16.*PI)
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
      NGAUS=MIN0((LPLUSX+10*(ICON(19)+1))/2,NPLM)
      CALL LEGAUS(2*NGAUS,XG,WG)
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
      C1=SQRT(T1-T2)/DQ(1)
      M1=MAX0(INT(C1),1)
      D1=C1-FLOAT(M1)
      C2=SQRT(T1+T2)/DQ(1)
      M2=MAX0(INT(C2),1)
      D2=C2-FLOAT(M2)
      C3=SQRT(T3-T4)/DQ(2)
      M3=MAX0(INT(C3),1)
      D3=C3-FLOAT(M3)
      C4=SQRT(T3+T4)/DQ(2)
      M4=MAX0(INT(C4),1)
      D4=C4-FLOAT(M4)
      MX=KRR*(M3-2)+M1
      MY=KRR*(M4-2)+M2
      DO 140 J=1,4
      FX(J)=0.0
      MX=MX+KRR
      FY(J)=0.0
      MY=MY+KRR
      IF(M1.lt.KRR-2 .and. M3.lt.KP -2 
     1               .and. M2.lt.KRR-2 .and. M4.lt.KP -2) then
C     INTERPOLATE IN Q DIRECTION
      FX(J)=
     1   -D1*(D1-1.)*(D1-2.)*FF(MX  )/6.+(D1**2-1.)*(D1-2.)*FF(MX+1)/2.  
     2   -D1*(D1+1.)*(D1-2.)*FF(MX+2)/2.+D1*(D1**2-1.)*FF(MX+3)/6.
C
      FY(J)=
     1   -D2*(D2-1.)*(D2-2.)*FF(MY  )/6.+(D2**2-1.)*(D2-2.)*FF(MY+1)/2.
     2   -D2*(D2+1.)*(D2-2.)*FF(MY+2)/2.+D2*(D2**2-1.)*FF(MY+3)/6.
C
      endif
  140 continue
C
C     INTERPOLATE IN P DIRECTION
      G1=
     1   -D3*(D3-1.)*(D3-2.)*FX(1)/6.+(D3**2-1.)*(D3-2.)*FX(2)/2.
     2   -D3*(D3+1.)*(D3-2.)*FX(3)/2.+D3*(D3**2-1.)*FX(4)/6.
C
      G2=
     1   -D4*(D4-1.)*(D4-2.)*FY(1)/6.+(D4**2-1.)*(D4-2.)*FY(2)/2.
     2   -D4*(D4+1.)*(D4-2.)*FY(3)/2.+D4*(D4**2-1.)*FY(4)/6.
C
  160 CONTINUE
      SUM1=SUM1+(G1+G2)*P(I,1)
      SUM3=SUM3+(G1+G2)*P(I,3)
      SUM2=SUM2+(G1-G2)*P(I,2)
      SUM4=SUM4+(G1-G2)*P(I,4)
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

      SUBROUTINE DSTRI2(IQ,DR,K,F,FR,QNUM,OPT,AR ,C,SI,SCALE)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      DIMENSION F(1025,2),FR( 4000),QNUM(4,2),TVCC(21),IQ(5)
     1,         D1(12),D2(12),AG(16),WG(16),AR(16),C(4),FLINT(24)
C
      Data PI/ 3.141592654/
c
      L =IQ(1)
      IS=IQ(2)
      JJ=IQ(3)
      IT=IQ(4)
      LP=IQ(5)
      LT=0
      DQ=C(1)*0.5
      KQ=C(2)
      KR=C(3)
      N1=QNUM(1,1)
      N2=QNUM(1,2)
      L1=QNUM(2,1)
      L2=QNUM(2,2)
      J1=QNUM(3,1)
      J2=QNUM(3,2)
      IS1=QNUM(4,1)
      IS2=QNUM(4,2)
      ITEMP=N1+N1+N2+N2
      NX=(ITEMP+L1+L2+L+2)/4+2
      NX=FLOAT((NX+5)/2)*(1.0+SCALE)
C     WRITE(6,9000)KR,NX
C9000 FORMAT(15H0 NO. R STEPS =,I3,18H     NO. X STEPS =,I3)
C
      NX=MIN0(NX,16)
      CALL LEGAUS(NX+NX,AG,WG)
C
C    CALCULATE D(LP) FOR LP
C
      IF(LP.EQ.0) GO TO 31
      TEMP=YXFCT(LP,2*LP)/SQRT(YXFCT(1,2*LP))
      DO 30 N=1,NX
      WG(N)=WG(N)*TEMP*(SQRT(1.0-AG(N)**2)/2.0)**LP
   30 CONTINUE
   31 CONTINUE
      LL=L+L
      LL1=L1+L1
      LL2=L2+L2
      LP2=LP+LP
      LT2=LT+LT
      FACT1=YXFCT(L1,LL1)/2.0**L1
      FACT2=YXFCT(L2,LL2)/2.0**L2
      FNORM=SQRT(FLOAT((2*LP+1)*(LL1+1)*(LL2+1))/FLOAT(LL+1))
      FNORM=FNORM/VCC(LL,LP2,LL,0,LP2)
      TEMP=2.0
      IZ=1
      DO 50 I=1,4
      IF(QNUM(I,1).NE.QNUM(I,2)) GO TO 55
   50 CONTINUE
      TEMP=4.0
      IZ=0
   55 CONTINUE
      FSYM=-PHASEF(IS/2+IT+LP)
      FNORM=FNORM*OPT/(SQRT(TEMP)*2.0*PI**2)
     1           *SQRT(FLOAT((LL+1)*(IS+1)*(J1+1)*(J2+1)))
     2           *WINEJ(LL1,IS1,J1,LL2,IS2,J2,LL,IS,JJ)
      IF(FNORM.EQ.0.0) GO TO 1000
C
C
      LPL2=MAX0(LL1,LL2)+2
      DO 70 M=2,LPL2
      FLINT(M)=1.0/FLOAT(M-1)
   70 CONTINUE
      FLINT(1)=0.0
      LPL2=MIN0(LL1,LL2,L1+L2-LP)+1
      M1=MIN0(L1,L2+LP)
      M2=LP-M1
      M1Z=M1
      M2Z=-M2
      DO 80 M=1,LPL2
      PHAS=1.0
      IF(M1.GT.0) PHAS=PHAS*PHASEF(M1)
      IF(M2.GT.0) PHAS=PHAS*PHASEF(M2)
      TVCC(M)=VCC(LL1,LL2,LL,2*M1,2*M2)*PHAS*PHASEF(M2)
      M1=M1-1
      M2=M2+1
   80 CONTINUE
      FM1=1.0
      FL1=L1
      LPL=L1+1
      DO 90 M=1,LPL
      M1=M1Z-M+2
      IF(M1.GT.0   ) TVCC(M1)=TVCC(M1)/SQRT(FM1)
      M1=M1Z+M
      IF(M1.LE.LPL2) TVCC(M1)=TVCC(M1)/SQRT(FM1)
      FM1=FM1*(FL1+FLOAT(M))*(FL1-FLOAT(M)+1.0)
   90 CONTINUE
      FM2=1.0
      FL2=L2
      LPL=L2+1
      DO 95 M=1,LPL
      M2=M2Z-M+2
      IF(M2.GT.0   ) TVCC(M2)=TVCC(M2)/SQRT(FM2)
      M2=M2Z+M
      IF(M2.LE.LPL2) TVCC(M2)=TVCC(M2)/SQRT(FM2)
      FM2=FM2*(FL2+FLOAT(M))*(FL2-FLOAT(M)+1.0)
   95 CONTINUE
C
C     CENTER OF MASS R LOOP
C
      INDEX=0
      R=0.0
      DO 500 M=1,KQ
      R=R+DQ
      QL=(2.0*R)**L
      RSQ=R**2
      SUMR=0.0
C
C     RELATIVE R LOOP
C
      DO 400 MR=1,KR
      INDEX=INDEX+1
      RS=AR(MR)
      RSSQ=RS**2+RSQ
      RPROD=2.0*R*RS
C
C     RELATIVE R ANGLE LOOP
C
      SUMX=0.0
      DO 300 MX=1,NX
      X=AG(MX)
      WX=WG(MX)
      IX=0
      TEMP=1.0
      R1=SQRT(RSSQ-RPROD*X)
      R2=SQRT(RSSQ+RPROD*X)
      FK1=R1/DR
      K1=MAX0(INT(FK1),1)
      FK1=FK1-FLOAT(K1)
      IF(K1.GT.K) GO TO 300
      FK2=R2/DR
      K2=MAX0(INT(FK2),1)
      FK2=FK2-FLOAT(K2)
      IF(K2.GT.K) GO TO 300
      COS1=(R-RS*X)/R1
      COS2=(R+RS*X)/R2
      SIN1=SQRT(1.0-COS1**2)
      SIN2=SQRT(1.0-COS2**2)
  120 CONTINUE
C
C    CALCULATE D(M) FOR L1
C
      DL1=FACT1
      IF(L1.EQ.0) GO TO 141
      DO 140 LI=1,L1
      DL1=DL1*SIN1
  140 CONTINUE
  141 CONTINUE
      D1(L1+1)=DL1
      D1(L1+2)=0.0
C
C     CALCULATE D(M) FOR L2
C
      DL2=FACT2
      IF(L2.EQ.0) GO TO 151
      DO 150 LI=1,L2
      DL2=DL2*SIN2
  150 CONTINUE
  151 CONTINUE
      D2(L2+1)=DL2
      D2(L2+2)=0.0
      IF(L1.EQ.0) GO TO 171
      FM1=LL1
      CTN1=COS1/SIN1
      DO 170 LI=1,L1
      D1(L1+1-LI)=(FM1*CTN1*D1(L1+2-LI)-D1(L1+3-LI))
     1            *FLINT(LI+1)*FLINT(LL1+2-LI)
      FM1=FM1-2.0
  170 CONTINUE
  171 CONTINUE
      IF(L2.EQ.0) GO TO 181
      FM2=LL2
      CTN2=COS2/SIN2
      DO 180 LI=1,L2
      D2(L2+1-LI)=(FM2*CTN2*D2(L2+2-LI)-D2(L2+3-LI))
     1            *FLINT(LI+1)*FLINT(LL2+2-LI)
      FM2=FM2-2.0
  180 CONTINUE
  181 CONTINUE
      PROD=0.0
      DO 185 LI=1,LPL2
      M1=IABS(M1Z-LI+1)+1
      M2=IABS(LI-1-M2Z)+1
      PROD=PROD+D1(M1)*D2(M2)*TVCC(LI)
  185 CONTINUE
  280 CONTINUE
C     INTERPOLATE - 4 POINT FORMULA
      F1=-FK1*(FK1-1.)*(FK1-2.)*F(K1-1,1)/6.
     1  +  (FK1**2-1.)*(FK1-2.)*F(K1  ,1)/2.
     2  - FK1*(FK1+1.)*(FK1-2.)*F(K1+1,1)/2.
     3  +       FK1*(FK1**2-1.)*F(K1+2,1)/6.
      F2=-FK2*(FK2-1.)*(FK2-2.)*F(K2-1,2)/6.
     1  +  (FK2**2-1.)*(FK2-2.)*F(K2  ,2)/2.
     2  - FK2*(FK2+1.)*(FK2-2.)*F(K2+1,2)/2.
     3  +       FK2*(FK2**2-1.)*F(K2+2,2)/6.
C     INTERPOLATE - 2 POINT FORMULA
C      F1=FK1*(F(K1+1,1)-F(K1,1))+F(K1,1)
C      F2=FK2*(F(K2+1,2)-F(K2,2))+F(K2,2)
      SUMX=SUMX+WX*PROD*F1*F2*TEMP
      IF(IX.NE.0) GO TO 300
      IX=1
      TEMP=FSYM
      ITEMP=K1
      K1=K2
      K2=ITEMP
      ATEMP=FK1
      FK1=FK2
      FK2=ATEMP
      IF(IZ.EQ.0) GO TO 280
      ATEMP=COS1
      COS1=COS2
      COS2=ATEMP
      ATEMP=SIN1
      SIN1=SIN2
      SIN2=ATEMP
      GO TO 120
  300 CONTINUE
      FR(INDEX)=FR(INDEX)+SUMX*FNORM/QL
      IF(M+MR.EQ.2) SI=SUMX*FNORM/QL
  400 CONTINUE
  500 CONTINUE
 1000 CONTINUE
C     Q=0.0
C     DO 1111 M=1,10
C     Q=Q+DQ*2.
C     WRITE(6,7778)Q,(FR(MR,M),MR=1,KR)
C7778 FORMAT(3H Q=,F5.2,1P10E12.5/(8X,1P10E12.5))
C     WRITE(6,7778)
C1111 CONTINUE
C     stop 'in bxdwk5'
      RETURN
      END
