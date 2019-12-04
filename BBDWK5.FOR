      SUBROUTINE BBDWK5
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      COMMON/BDWCK_5/LTR,LTP,JTR,JTP,T(4),DQ(2),QMAX(2)
c
      CALL RADINT(LTR,LTP,JTR,JTP,T)
      RETURN
      END


      SUBROUTINE RADINT(LTR,LTP,JTR,JTP,T)
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
      dimension nadd(2*nls), ak1(imat,16), bk1(imat,16), BICO(16,2)
     1,buffer(ispc1), store(imat*imat), fll(ispc1/2)
     2,T(4), SUM(24), C(2)
c
      equivalence (space0(   1),nadd(1)),  (space0( 401),ak1(1,1))
     1,           (space0(1401),bk1(1,1)), (space0(2501),store(1))
     2,           (space1(   1),fll(1)),   (space2(   1),buffer(1))
C
      NBUFF=IFLL*2
      LTR2=LTR+LTR
      LTP2=LTP+LTP
      JR=NS(1)
      JS=NS(2)
C     I1=1 ==> Pickup    I1=2 ==> Stripping
      if    (fma(1) .gt. fma(2)) then
        I1=1
      else
        I1=2
      endif
C     I2=1 ==> Post  stripping, prior pick up 
C     I2=2 ==> Prior stripping, post  pick up 
      fAmass = min(fma(1), fma(2))
      fbmass = min(fm (1), fm (2)) 
      if    (abs(fAmass - fbmass) .lt. 0.20) then
c     special case if projectile core and target core masses are the same 
      if(i1 .eq. 1 )then 
        if    (ecm(3) .gt. ecm(4)) then
          I2=2
        else
          I2=1
        endif
      else
        if    (ecm(3) .gt. ecm(4)) then
          I2=1
        else
          I2=2
        endif 
      endif
      elseif(abs(fma(3) - fAmass) .lt. 0.20) then
c     post  stripping, prior pick up 
        I2=1
      elseif(abs(fma(3) - fbmass) .lt. 0.20) then 
c     prior stripping, post  pick up 
        I2=2
      else 
        go to 4000
      endif
c
      LTMIN=MAX0(IABS(LTR-LTP),IABS(JTR-JTP)/2)+1
      LTMAX=MIN0(     LTR+LTP ,    (JTR+JTP)/2)+1
      LPLUSX=LPLUS+LTMAX-1
C
C
C      CALL KSECOND(TIME)
C      WRITE(6,9950)TIME
C
C     CONSTRUCT MOMENTUM TABLES
C
      LSQ=MAX0(LTR,LTP)+1
      DO 510 N=1,NMAX
      IF(EZER( N).EQ.0.0) GO TO 510
      AK1(N,1)=1.0
      IF(LSQ .ne. 1) then
      DO 509 I=2,LSQ
      AK1(N,I)=AK1(N,I-1)*EZER(N)
  509 CONTINUE
      endif
  510 CONTINUE
      DO 513 N=1,NMAX
      NN=N+NMAX
      IF(EZER(NN).EQ.0.0) GO TO 513
      BK1(N,1)=1.0
      IF(LSQ .ne. 1) then
      DO 512 I=2,LSQ
      BK1(N,I)=BK1(N,I-1)*EZER(NN)
  512 CONTINUE
      endif
  513 CONTINUE
      if(i2.eq.1) then
      IF(I1.EQ.1) then
C     PICKUP
        LX=LTP+1
        LY=LTR+1
        TX=T(4)
        TY=T(1)
      else
C     STRIPPING
        LX=LTR+1
        LY=LTP+1
        TX=T(2)
        TY=T(3)
      endif
      else
      IF(I1.EQ.1) then
C     PICKUP
        LX=LTR+1
        LY=LTP+1
        TX=T(2)
        TY=T(3)
      else
C     STRIPPING
        LX=LTP+1
        LY=LTR+1
        TX=T(4)
        TY=T(1)
      endif
      endif
	write(6,'(a,4i5,2f10.4)')' I1, I2, Lx, Ly, TX, TY'
     1                           , I1, I2, Lx, Ly, TX, TY
c
      LTX2=LX+LX-2
      LTY2=LY+LY-2
C     INITIAL STATE
      BICO(1,1)=1.0
      IF(LX .ne. 1) then
      T1=LX+LX-2
      DO 517 I=2,LX
      BICO(I,1)=-BICO(I-1,1)*SQRT(T1*(T1-1.0)/FLOAT((2*I-3)*(2*I-2)))*TX
      T1=T1-2.0
  517 CONTINUE
      endif
C     FINAL STATE
      BICO(1,2)=1.0
      IF(LY .ne. 1) then
      T1=LY+LY-2
      DO 519 I=2,LY
      BICO(I,2)=-BICO(I-1,2)*SQRT(T1*(T1-1.0)/FLOAT((2*I-3)*(2*I-2)))*TY
      T1=T1-2.0
  519 CONTINUE
      endif
      LTOT=LTR+LTP
      NLTX=2*LTOT+1
C
C     CHECK ON DW BUFFER STORAGE CAPACITY
C
      I=0
      DO 525 L=1,NLTX
      I=I + ((N2(L,1)-N1(L,1)+1)*NS(1) + (N2(L,2)-N1(L,2)+1)*NS(2))*2
  525 CONTINUE
      IF(I.GT.NBUFF) THEN
      WRITE(6,9102)I,NBUFF
 9102 FORMAT('0DISTORTED WAVE BUFFER TOO SMALL, NEED ',I5
     1,' , ',I5, ' AVAILABLE')
      stop 'Distorted wave buffer too small'
      ENDIF
C
C     L TRANSFER LOOP
C
      DO 2500 LT=LTMIN,LTMAX
      LZ=LT-1
      LX2=LZ+LZ
C
C      CALL KSECOND(TIME)
C      WRITE(6,9951)TIME,LZ
C
      N=LTR+LTP-LZ
      IPHAS=N-2*(N/2)
      NLTX=LTOT+LTOT+1
      INC=LX2+2-IPHAS-IPHAS
      INCR=INC*LPLUS
      JRS=NS(1)*NS(2)
      IFL=INCR*JRS
      if(ifl .gt. ispc1/2) write(6,'(a,i6,a,i6)')
     1' Radial matrix element storage exceeds',ispc1/2,' --> needs',ifl
      JRS=JRS+JRS
C
C     CLEAR STORAGE FOR RADIAL INTEGRALS
C
      DO 535 I=1,IFL
      FLL(I)=0.0
  535 CONTINUE
      LL=0
      NA=0
      NB=0
C
C     FORM FACTOR EXPANSION LOOP
C
      DO 2000 LK=1,LPLUSX
      L=LK-1
      LL2=LK+LK-2
      FACTL=LL2+1
C
C     FILL BUFFER WITH DISTORTED WAVES
C
      IF(L.EQ.0   ) GO TO 540
      IF(L.LE.LTOT) GO TO 550
  540 CONTINUE
      LL=LL+1
      IF(LL.GT.LPLUS) GO TO 550
  542 CONTINUE
      NADD(LL    )=NA
      NB=NA+(N2(LL,1)-N1(LL,1)+1)*NS(1)*2
      NADD(LL+200)=NB
      NA=NB+(N2(LL,2)-N1(LL,2)+1)*NS(2)*2
      IF(NA.LT.NBUFF) GO TO 545
      IF(NLTX.GT.0) GO TO 544
      NA=0
      GO TO 542
  544 CONTINUE
      WRITE(6,9101)LZ
 9101 FORMAT(49H0DISTORTED WAVE BUFFER TOO SMALL FOR L TRANSFER =,I3)
      stop 'Distorted wave buffer too small'
  545 CONTINUE
      J1=NADD(LL    )+1
      J2=NB
      J3=NB+1
      J4=NA
      READ (4)(BUFFER(J),J=J1,J2)
      READ (4)(BUFFER(J),J=J3,J4)
      NLTX=NLTX-1
      IF(NLTX.GT.0) GO TO 540
  550 CONTINUE
C
C     READ IN LEGENDRE EXPANSION COEFFICENTS
      READ (3)INDEX,NIA,NIB,NFA,NFB,(STORE(I),I=1,INDEX)
C
C
      LP1=MAX0(L-LTOT,0)+1
      LP2=MIN0(L+LTOT+1,LPLUS)
C
C     INITIAL WAVE L LOOP
C
      DO 1600 L1=LP1,LP2
      LI=L1-1
      LI2=LI+LI
      NI1=N1(L1,1)
      NI2=N2(L1,1)
      LTEMP=LTOT-IABS(L-LI)
      LP3=IABS(LI-LZ)+1+IPHAS
      LP4=     LI+LZ +1-IPHAS
C
C     FINAL   WAVE L LOOP
C
      DO 1500 L2=LP3,LP4,2
      LF=L2-1
      IF(IABS(LF-L).GT.LTEMP) GO TO 1500
      IF(L2.LT.LP1) GO TO 1500
      IF(L2.GT.LP2) GO TO 1500
      LF2=LF+LF
      IADD=LF*INC+LI-LF+LZ-IPHAS+1
      NF1=N1(L2,2)
      NF2=N2(L2,2)
C
      DO 575 I=1,JRS
      SUM(I)=0.0
  575 CONTINUE
      LAMIA=IABS(L-LI)+1
      LAMIB=MIN0(L+LI,LTR+LTP)+1
C
C     LAMI LOOP
C
      DO 1200 LAMI=LAMIA,LAMIB,2
      LAMI2=LAMI+LAMI-2
      TEMP1=VCC(LAMI2,LL2,LI2,0,0)*SQRT(FLOAT(LAMI2+1))*PHASEF(LAMI-1)
      LAMFA=IABS(L-LF)+1
      LAMFB=MIN0(L+LF,LTR+LTP)+1
C
C     LAMF LOOP
C
      DO 1100 LAMF=LAMFA,LAMFB,2
      I1=LAMI+LAMF-2
      IF(I1.LT.LZ  ) GO TO 1100
      IF(I1.GT.LTOT) GO TO 1100
      LAMF2=LAMF+LAMF-2
      TEMP2=TEMP1*VCC(LAMF2,LL2,LF2,0,0)*SQRT(FLOAT(LAMF2+1))
     1*RACAH(LI2,LAMI2,LF2,LAMF2,LL2,LX2)
C
      DO 900 LMX=1,LX
      LMX2=LMX+LMX-2
      JK=NMAX*(LMX-1)+NF1-1
      JX=NMAX*(LX-LMX)+NI1-1
      TEMP3=TEMP2*BICO(LMX,1)
C
      DO 800 LMY=1,LY
      LMY2=LMY+LMY-2
      IF(PHASEF(LY-LMY+LMX+LAMF  ).LT.0.0) GO TO 795
      IF(PHASEF(LX-LMX+LMY+LAMI  ).LT.0.0) GO TO 795
      IF(LTY2-LMY2+LMX2.LT.LAMF2) GO TO 795
      IF(LTX2-LMX2+LMY2.LT.LAMI2) GO TO 795
      JL=NMAX*(LMY-1)+NI1-1
      JY=NMAX*(LY-LMY)+NF1-1
      TEMP =TEMP3*BICO(LMY,2)
     1*VCC(LTY2-LMY2,LMX2,LAMF2,0,0)*VCC(LTX2-LMX2,LMY2,LAMI2,0,0)
     2*WINEJ(LX2,LTY2,LTX2,LAMF2,LTY2-LMY2,LMX2,LAMI2,LMY2,LTX2-LMX2)
      INDA=NFB-NFA+1
C
      IK=0
      DO 750 ISI=1,JR
      DO 749 ISF=1,JS
      ISTOR=INDA*(NI1-NIA)-NFA+1
      IADD1=NADD(L1    )+ISI+ISI-1
      IND2=JL
      IND3=JX
      SUM1=0.0
      SUM2=0.0
C     SUM OVER INITIAL STATE MOMENTA
      DO 700 NI=NI1,NI2
      IADD2=NADD(L2+200)+ISF+ISF-1
      IND2=IND2+1
      IND3=IND3+1
      FACT1=TEMP*AK1(IND2,1)*AK1(IND3,1)
      C(1)=FACT1*BUFFER(IADD1  )
      C(2)=FACT1*BUFFER(IADD1+1)
      IND1=JK
      IND4=JY
C     SUM OVER FINAL STATE MOMENTA
      DO 600 NF=NF1,NF2
      IND1=IND1+1
      IND4=IND4+1
      ELEM=BK1(IND1,1)*BK1(IND4,1)*STORE(ISTOR+NF)
      SUM1=SUM1
     1    +(C(1)*BUFFER(IADD2  )-C(2)*BUFFER(IADD2+1))*ELEM
      SUM2=SUM2
     1    +(C(2)*BUFFER(IADD2  )+C(1)*BUFFER(IADD2+1))*ELEM
      IADD2=IADD2+JS+JS
  600 CONTINUE
      IADD1=IADD1+JR+JR
      ISTOR=ISTOR+INDA
  700 CONTINUE
      IK=IK+1
      SUM(IK+IK-1)=SUM(IK+IK-1)+SUM1*FACTL
      SUM(IK+IK  )=SUM(IK+IK  )+SUM2*FACTL
  749 CONTINUE
  750 CONTINUE
  795 CONTINUE
  800 CONTINUE
  900 CONTINUE
 1100 CONTINUE
 1200 CONTINUE
      INDEX=IADD
      IK=0
      DO 1450 ISI=1,JR
      DO 1449 ISF=1,JS
      IK=IK+1
      FLL(INDEX  )=FLL(INDEX  )+SUM(IK+IK-1)
      FLL(INDEX+1)=FLL(INDEX+1)+SUM(IK+IK  )
      INDEX=INDEX+INCR
 1449 CONTINUE
 1450 CONTINUE
 1500 CONTINUE
 1600 CONTINUE
 2000 CONTINUE
      REWIND 3
      REWIND 4
C
C     WRITE RADIAL MATRIX ELEMENTS ON TAPE 2
C
      IND=0
      DO 2120 I=1,JR
      DO 2110 J=1,JS
      INC=IND+1
      IND=INC-1+INCR
      WRITE(2)(FLL(II),II=INC,IND)
 2110 CONTINUE
 2120 CONTINUE
      IF(ICON(7).ne.0) then
C
C     WRITE( RADIAL MATRIX ELEMENTS IF ICON(7).NE.0
C
      IND=0
      IS1=-IS(1)
      DO 2200 I=1,JR
      IS2=-IS(2)
      DO 2198 J=1,JS
      WRITE(6,9999)ALPHA,(IDAT(ID),ID=1,3)
      WRITE(6,9900)IS2,IS1
      WRITE(6,9901)
      DO 2190 LL=1,LPLUS
      LM=LL-1
      INC=IND+1
      IND=IND+LX2+2-IPHAS-IPHAS
      WRITE(6,9902)LM,(FLL(INDEX),INDEX=INC,IND)
 2190 CONTINUE
      IS2=IS2+2
 2198 CONTINUE
      IS1=IS1+2
 2200 CONTINUE
      endif
      TEMP=RACAH(2*LTR,2*LTP,JTR,JTP,LX2,IS(3))
      T1=DSIGN(1.0D0,TEMP)
      TEMP=FLOAT(LX2+1)*TEMP**2
      WRITE(6,9100)LTR,LTP,JTR,JTP,LZ,IS(3),TEMP,T1
 2500 CONTINUE
      RETURN
c
 4000 continue
      write(6,9910)
      ibf(3) = 1
      return 
c
 9900 FORMAT(24H0 RADIAL MATRIX ELEMENTS ,9H,  J2=L2+,I2,2H/2,9H,  J1=L1
     1+,I2,2H/2 )
 9100 FORMAT('0(2*P+1)*RACAH(LTR,LTP,JTR,JTP,P,IS)**2 FOR '
     1,2I3,2(I3,2H/2),2I3,4H/2 =,F8.5,'   SIGN =',F3.0)
 9901 FORMAT('0 L2  F(L2, L2-LTR   )    F(L2, L2-LTR +2)    '
     1,            'F(L2, L2-LTR +4)    F(L2, L2-LTR +6)')
 9902 FORMAT(I4,1P,10E10.3/(4X,1P,10E10.3))
 9910 format('0 Error, in input, inconsistent masses in RADINT ')
 9950 FORMAT('0ELAPSED TIME =',F10.3,'     SUBR. RADINT BEGUN')
 9951 FORMAT('0ELAPSED TIME =',F10.3,'     R.M.E. BEGUN FOR LTR =',I2)
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END
