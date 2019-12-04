c$debug
      SUBROUTINE ADWCK5
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical i_sym(2)
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
      dimension dtemp(2400), plm(1600), sigplt(200,2), d(2000), rsig(2)
c
      equivalence (Space0(   1),dtemp(1)),  (Space0(2401),plm(1))
     1,           (Space0(4001),sigplt(1,1))
     2,           (Space1(   1),d(1))
c
      CALL AADWK5
c
c     check for identical particles      
c
      do 40 i=1,2
      if(fm(i).eq.fma(i).and.z(i).eq.za(i)) then
        i_sym(i)=.true.
      else
        i_sym(i)=.false.
      endif
   40 continue 
c
      CALL INTEG5(i_sym,rsig)
c
      if(icon(6).ne.0.or.icon(16).ne.0) then
C
C     PRINT ELASTIC CROSS-SECTIONS IF ICON(6).NE.0
C
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      WRITE(6,9904)
C
C
C     CALL OVERLAY(FNAME,1,3)
        CALL ELSIG(dtemp,d,plm,sigplt,angle,fk,eta,rsig,alpha
     1            ,idat,is,icon,lplus,i_sym)
C
      endif
c     Terminate program here if icon(20).eq.1
      if(icon(20).eq.1) then
        IBF(3)=1
        WRITE(6,9300)
      endif
      RETURN
c
 9300 FORMAT('0PROGRAM TERMINATED BY ICON(20)')
 9904 FORMAT('0ELASTIC SCATTERING CROSS-SECTION')
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

      SUBROUTINE AADWK5
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
c
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (eta2 = 0.10, eta7 = 1.0e-12)
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
      common/tnt/ pq(4,8),fmtemp(16)
c
      dimension u(1600),v(1600), ub(1600), space(1600)
     1,         f(nl2),fp(nl2),g(nl2),gp(nl2),s(nl2)
     2,         scale(8),q(4),p(4),title(10), ff(10) 
c
      equivalence (Space0(   1),u(1)), (Space0(1601),v(1))
     1,           (Space0(3201      ),f(1)), (Space0(3201+  nl2),fp(1))
     2,           (Space0(3201+2*nl2),g(1)), (Space0(3201+3*nl2),gp(1))
     3,           (Space0(3201+4*nl2),s(1))
     4,           (Space1(   1),ub(1)), (Space1(1601),space(1))
c
      LAM=0
      IF(ICON(1).EQ.0) GO TO 5
C
C     READ CARD SET 2
      READ (5,9001)ANGLE
      IF(ICON(1).EQ.1) GO TO 5
      II=(ANGLE(1)-ANGLE(2))/ANGLE(3)+1.0
      IF(II.LT.0) II=0
      ANGLE(1)=II
    5 CONTINUE
      WRITE(6,9010)ANGLE
C
C
C     READ CARD SET 3
      READ (5,9002)Linput,NLTR,(SCALE(I),I=1,NLTR)
      DO 8  II=1,NLTR
      IF(SCALE(II).EQ.0.0) SCALE(II)=1.0
    8 CONTINUE
      WRITE(6,9008)Linput,NLTR,(SCALE(I),I=1,NLTR)
      MLTR=NLTR
C
C
C     READ CARD SET 4
      READ (5,9001)DRF,RMAX,ACC
      WRITE(6,9011)DRF,RMAX,(ACC(I),I=1,6)
C
      IF(DRF.EQ.0.0) DRF=0.1
      K =ABS(RMAX)/DRF+1.0E-08
      K=MIN0(K,KMAX)
      KC=2*((K+1)/2)
      LPLUS=IABS(Linput)+1
      LC=LPLUS
      IBF(5)=0
      IBF(6)=0
      WRITE(6,9504)
C
C     READ IN DISTORTED WAVE INFORMATION
C
      do 10 i=1,2
c     ibf(4) is spin orbit potential flag,       
      IBF(4)=0
      ind=(i-1)*2*maxk + 1
      CALL FORMF(U(ind),V(ind),i,0,Q, ff)
      if(IBF(4).NE.0) then
        IBF(i+4)=1
      else  
        IBF(i+4)=0
      endif
      
C
C       PRINT OUT K(R)**2 FOR DISTORTED WAVES
C 
      if(icon(15) .ne. 0) then
        write(6,9057)i
        call ffprint(kc,icon(15),u(ind),dr(i))
      if(ibf(i+4) .ne. 0)then
        write(6,9058)i
        call ffprint(kc,icon(15),v(ind),dr(i))
      endif
      endif
c      
   10 continue
c        
      IND=0
      INF=0
C
C     FORM FACTOR LOOP
C
      DO 200 II=1,NLTR
      ICC = 1
      WRITE(6,9999)ICC,ALPHA,(IDAT(I),I=1,3)
      WRITE(6,9503)
C     IF(SCALE(II).EQ.0.0) SCALE(II)=1.0
      IF(ICON(2).NE.0) GO TO 100
C
C     FIRST FORM FACTOR
C
c
      do 20 i=1,10
      ff(i) = 0.0
   20 continue 
      IF(RMAX.LT.0.0) K=KC
      CALL FORMF(UB,SPACE(1)   ,3,LAM,Q, ff)
      LTRT(II)=Q(2)
      JTRT(II)=Q(3)
      ISTRT(II)=Q(4)
      KFF=KZ
C
      KT=KFF
      if(icon(4).ne.0) then
        call ffprint(kt,icon(4),ub,dr(3))
      endif
C
      KT2 = 2*KT
      WRITE(9)KT2,Q(1),ETA(3),FK(3),(UB(I),I=1,KT2), ff 
C
      LTRP(II)=0
      JTRP(II)=0
      ISTRP(II)=0
C     TEST FOR ZERO RANGE AND READ IN OPTIONS
      IF(ICON(11).NE.0) GO TO 35
      if(icon(4) .ne. 0) ICC = 0
      WRITE(6,9999)ICC,ALPHA,(IDAT(I),I=1,3)
      WRITE(6,9503)
      IBF(2)=0.0
C
C     SECOND FORM FACTOR
C
c
      do 30 i=1,10
      ff(i) = 0.0
   30 continue 
      CALL FORMF(UB,SPACE(1)   ,4,LAM,Q, ff)
      LTRP(II)=Q(2)
      JTRP(II)=Q(3)
      ISTRP(II)=Q(4)
      KFF=KZ
      IF(FK2(4).NE.0.0) THEN
C
C     IBF(2) = # points for Rmatch
C     calculate # points to Rmatch + 2fm.
C
      IM=2.0001/DRF
      IM=IM+IM+IBF(2)+IBF(2)-1
      FACT=SCALE(II)*HBARC**2/(2.0*AMU*FMU(4))
      UMAX=0.0
c
      DO 33 M=1,KFF
      MK=M+M-1
      UB(MK  )=UB(MK  )*(FK2(4)-SPACE(M))*FACT
      IF(FK2(4).LT.0.0) THEN
        UMAX= MAX (UMAX, ABS(UB(MK  )))
        IF(MK.GT.IM.AND. ABS(UB(MK  ))/UMAX.LT.ETA7) then
          KT = M
          go to 34
        ENDIF
      ENDIF
   33 CONTINUE
      ENDIF
   34 continue
      KT2 = 2*KT
      WRITE(9)KT2,Q(1),ETA(4),FK(4),(UB(I),I=1,KT2), ff 
      if(icon(4).eq.0) then
        WRITE(6,9034)
        call ffprint(kt,icon(4),ub,dr(4))
      endif
      GO TO 160
C
C     HERE FOR ZERO RANGE OR READ IN OPTIONS FOR FOURIER TRANSFORM
C
   35 CONTINUE
      READ (5,9041)Q,DQ,TITLE
      IF(DQ.EQ.0.0) DQ=0.05
      IMAX=    Q(1)
      LTRP(II)=Q(2)
      JTRP(II)=Q(3)
      ISTRP(II)=Q(4)
      WRITE(6,9507)TITLE
      WRITE(6,9508)IMAX,LTRP(II),JTRP(II),ISTRP(II),DQ
      IF(ICON(11).EQ.1) THEN
C     ZERO RANGE FORM FACTOR
C     RQ IS FINITE RANGE PARAMETER
      RQ=Q(1)
      WRITE(6,9007)RQ
      IMAX=256
      QQ=0.0
      DO 37 I=1,IMAX
      QQ=QQ+DQ
      UB(I)=-100.*SCALE(II)/(1.0+(QQ*RQ)**2)
   37 CONTINUE
      ELSE
C     READ IN FORM FACTOR
      READ (5,9003)(UB(M),M=1,IMAX)
      DO 45 M=1,IMAX
      UB(M) =UB(M) *SCALE(II)
   45 CONTINUE
      ENDIF
      WRITE(9)IMAX,DQ,ETA(4),FK(4),(UB(I),I=1,IMAX)
      GO TO 160
C
C     TNT FORM FACTOR
C
C     ICON(11).EQ.0  Default   COMPUTE TARGET FF, READ IN PROJ.FF
C     ICON(11).EQ.1  ZR OPTION COMPUTE TARGET FF, ZR PROJ.FF a'la DWUCK4
C     ICON(11).EQ.2  READ IN   (TARGET*PROJ.) FF  FROM TAPE8 (see bxdwk5)
C     ICON(11).EQ.3  READ IN    TARGET FF
C     ICON(11).EQ.4  SAVE (TARGET*PROJ.) FF ON TAPE 8 FOR '2'
C     ICON(11).EQ.5  SAVE (TARGET*PROJ.) FF ON TAPE 8 FOR '3'
C
  100 CONTINUE
      ICNT=1
      IF(ICON(11).gt.5) then
        stop 'Invalid value for icon(11)'
      endif
c
      IF(ICON(11).EQ.2) GO TO 158
c
      READ (5,9000)Q,T1
      LTRT(II)=Q(2)
      JTRT(II)=Q(3)
      ISTRT(II)=Q(4)
      NT1=T1
      IF(NT1.NE.0) THEN
      DO 108 I=1,NT1
      READ (5,9050)TITLE
      WRITE(6,9050)TITLE
  108 CONTINUE
      ENDIF
      IF(ICON(11).EQ.3) GO TO 121
C     READ IN TARGET FORM FACTOR INFORMATION
  110 CONTINUE
C
C     CNTRL  = NO. OF S.P. ORBITALS
C     QCODE  = NOT USED ( FOR DWUCK4 COMPATIBILITY)
C     FMUV   = RMS RADIUS OF TRITON
C     OPT    = CONFIGURATION AMPLITUDE
C     TSPIN  = ISOSPIN OF TRANSFERED PARTICLES
C
      READ (5,9000)CNTRL,QCODE,FMUV,OPT,TSPIN
      WRITE(6,9006)ICNT,CNTRL,OPT,TSPIN
      ICN=ABS(CNTRL)
      ICNT=ICNT+ICN
      DO 120 IC=1,2
      IF(ICN.NE.0) GO TO 114
      KFF=1
      GO TO 116
  114 CONTINUE
      IF(ICN.EQ.1.AND.IC.EQ.2) GO TO 120
      CALL FORMF(SPACE(1),SPACE(0801),3,0,Q, ff)
      KFF=KZ
  116 CONTINUE
      WRITE(9) Q,Z(3),TSPIN,CNTRL,OPT,DRF,KFF,(SPACE(2*I-1),I=1,KFF)
      IF(ICN.EQ.0) GO TO 121
  120 CONTINUE
  121 CONTINUE
      IF(CNTRL.GT.0.0) GO TO 110
C
C     READ IN PROJECTILE FORM FACTOR
      READ (5,9000)Q,DEL
C
      IF(ICON(11).EQ.1.OR.ICON(11).EQ.5) GO TO 123
      GO TO 141
  123 CONTINUE
C
C     ZERO RANGE OPTION FOR TWO NUCLEON TRANSFER
C     GAUSSIAN FORM - SAME AS IN DWUCK4
C     RP     = FINITE RANGE PARAMETER
      RP=Q(1)
      WRITE(6,9007)RP
C
      IF(DEL.EQ.0.0) DEL=1.7
      X=10.0
      DX=0.100
      P(1)=DX/DEL
      P(3)=AINT(X/DX+1.E-6)
      P(2)=0.2
      P(4)=40.0
      WRITE(9)P
      KR=P(3)
      KP=P(4)
c     Fnorm contains the Dwuck4 "9.72" normalization.
      FNORM=-100.0*SCALE(II)*(8.0*3.14159*DEL**2)**0.75
c     nn c.m. momentum loop
      PP=0.0
      DO 140 M=1,KP
      PP=PP+P(2)
      T1=FNORM/(1.0+(PP*RP)**2)
c     nn rel. momentum loop
      QQ=0.0
      DO 130 N=1,KR
      QQ=QQ+P(1)
      UB(N)=T1   *EXP(-(DEL*QQ)**2)
  130 CONTINUE
      WRITE(9)(UB(N),N=1,KR)
  140 CONTINUE
      GO TO 151
C
C     READ IN PROJECTILE FORM FACTOR
C
  141 CONTINUE
      WRITE(6,9042)
      READ (5,9009)P,TITLE
      WRITE(6,9019)P,TITLE
      WRITE(9)P
      KR=P(3)
      KP=P(4)
C
      DO 150 M=1,KP
      READ (5,9003)(UB(N),N=1,KR)
      DO 148 N=1,KR
      UB(N)=UB(N)*SCALE(II)
  148 CONTINUE
      WRITE(9)     (UB(N),N=1,KR)
  150 CONTINUE
  151 CONTINUE
c     nn c.m. momentum
      PQ(1,II)=P(2)
      PQ(2,II)=P(4)
c     nn rel. momentum
      PQ(3,II)=P(1)
      PQ(4,II)=P(3)
C
      LTRP(II)=Q(2)
      JTRP(II)=Q(3)
      ISTRP(II)=Q(4)
      if(icon(11).eq.4.or.icon(11).eq.5) then
      WRITE(8)LTRP(II),JTRP(II),ISTRP(II)
     1       ,LTRT(II),JTRT(II),ISTRT(II)
      endif
c     remainder of form factor written out in bxdwk5
      GO TO 159
  158 CONTINUE
C
C     read in (target*projectile) form factor
C
      READ (8)LTRP(II),JTRP(II),ISTRP(II)
     1       ,LTRT(II),JTRT(II),ISTRT(II)
c     remainder of form factor read in bxdwk5
c
  159 CONTINUE
      FMA(3)= MIN (FMA(1),FMA(2))
      FMA(4)= MIN (FM (1),FM (2))
      FM(3)=2.0
      FM(4)=2.0
C
C     END OF PROJECTILE FORM FACTORS
C
  160 CONTINUE
      IND=MAX0(IND,MIN0(LTRT(II)+LTRP(II),(JTRT(II)+JTRP(II))/2))
      INF=MAX0(INF,JTRT(II)/2)
C
C     CHECK ON SPINS AND STATISTICS
C
      IF(ISTRT(II).EQ.ISTRP(II)) GO TO 170
      IBF(3)=1
      WRITE(6,9907)ISTRT(II),ISTRP(II)
  170 CONTINUE
C
C     STORAGE CHECK ON MAX. LTR IN RADINT
C
      I1=MAX0(LTRT(II),LTRP(II))
      IF(I1 .ge. ICO) then
      IBF(3)=1
      WRITE(6,9908)ICO-1
      endif
      IF(PHASEF( IS(1)+IS(2)+ISTRT(II) ).GT.0.0) GO TO 180
      IBF(3)=1
      WRITE(6,9906)IS(1),IS(2),ISTRT(II)
  180 CONTINUE
      FMTEMP(2*II-1)=FMA(3)
      FMTEMP(2*II  )=FMA(4)
c     end of form factor loop
  200 CONTINUE
C
      IF(NLTR.GT.8) IBF(3)=1
      l=iabs(linput)
      WRITE(6,9028)L     ,(LTRT(I) ,I=1,NLTR)
      WRITE(6,9029)NLTR  ,(JTRT(I) ,I=1,NLTR)
      WRITE(6,9033)       (ISTRT(I),I=1,NLTR)
      WRITE(6,9030)       (LTRP(I) ,I=1,NLTR)
      WRITE(6,9031)       (JTRP(I) ,I=1,NLTR)
      WRITE(6,9032)       (ISTRP(I),I=1,NLTR)
C
C     Shut off spins if no spin orbit
C     ibf(5) and ibf(6) are spin orbit flags for distorted waves
c     
      is_off=IBF(5)+IBF(6)
      IBF(7)=IS(1)
      IBF(8)=IS(2)
      IF(is_off.eq.0) then
        IS(1)=0
        FS(1)=0.0
        NS(1)=1
        IS(2)=0
        FS(2)=0.0
        NS(2)=1
      endif
  250 CONTINUE
C
C     CHECK ON COMPATIBILITY OF LMAX AND STORAGE
C
      IF(Linput.LT.0) LPLUS=LC
C
      INF=2*(INF+1)*NS(1)*NS(2)
      INE=2*(IND+1)*NS(1)*NS(2)
      I1=11400
      LPLUS=MIN0(LPLUS,LMAX)
      L1=MIN0(I1/(INF+2*(IND+1)),(I1-6100)/INE)
      IF(LPLUS.GT.L1) WRITE(6,9004)
      IF(RMAX.LT.0.0) K=KC
      K=5*(K/5)
      RMAX=FLOAT(K)*DRF
      L2=INT( MIN (FK(1),FK(2))*RMAX)+10
      IF(LPLUS.GT.L2) WRITE(6,9005)
      LPLUS=MIN0(LPLUS,L1,L2)
C
      
      L=LPLUS-1
      LPL2=LPLUS+LPLUS
      WRITE(6,9506)L,K,NLTR
      REWIND 9
C
C
      IY=1
      DO 400 N=1,2
      X2=DRF*FK(N)*FLOAT(K)
      X1=X2-DRF*FK(N)
      CALL COU(X1,X2,ETA(N),LPLUS,ETA2,F(IY),FP(IY),G(IY),GP(IY),S(IY))
      IY=IY+LPLUS
  400 CONTINUE
C
C
      RETURN
C
 9000 FORMAT(10F8.4)
 9001 FORMAT(10F8.4)
 9002 FORMAT(2I3,2X,8F8.4)
 9003 FORMAT(5E16.7)
 9004 FORMAT(8H0WARNING,8(1H*),2X,30HLMAX DECREASED TO FIT STORAGE   )
 9005 FORMAT(8H0WARNING,8(1H*),2X,30HLMAX DECREASED TO FIT MAX K*R   )
 9006 FORMAT('0NO.',I4,'  CNTRL=',F6.2,6H  OPT=,F8.4,8H  TSPIN=,F6.2)
 9007 FORMAT('0FINITE RANGE PARAMETER =',F6.2)
 9008 FORMAT(18H0CARD SET 3 DATA  ,9H   LMAX =,I4,14H         NLTR=,I4
     1,/18X,18H    SCALE FACTORS=,8F9.4)
 9009 FORMAT(4F8.4,8X,10A4)
 9010 FORMAT(18H0ANGLE DATA        ,9H   THETN=,F9.4,9H   THET1=,F9.4
     1,9H   DTHET=,F9.4,9H   A-ANG=,F9.4,9H   B-ANG=,F9.4)
 9011 FORMAT(18H0CARD SET 4 DATA   ,9H   DRF  =,F9.4,9H   RMAX =,F9.4
     1,20H  ACCURACY SETTINGS, ,6F9.3)
 9019 FORMAT(8H0   DQ =,F6.2,4X,4HDP =,F6.2,4X,4HNQ =,F6.2
     1,4X,4HNP =,F6.2,8X,10A4)
 9028 FORMAT(18H0TARGET TRANSFERS ,9H   LMAX =,I4,8X,6HLTRT =,8I4)
 9029 FORMAT(21X,6HNLTR =,I4,7X,7H2*JTRT=,8I4)
 9030 FORMAT('0PROJECTILE TRANSFERS ',17X,6HLTRP =,8I4)
 9031 FORMAT(38X,7H2*JTRP=,8I4)
 9032 FORMAT(38X,7H2*ISTP=,8I4)
 9033 FORMAT(38X,7H2*ISTRT,8I4)
 9034 FORMAT('0FORM FACTOR * V(R)',/)
 9041 FORMAT(5F8.4,10A4)
 9042 FORMAT('0PROJECTILE FORM FACTOR DATA')
 9050 FORMAT(10A4)
 9057 FORMAT('0Central    potential K(R)**2 FOR DISTORTED WAVE',I2)
 9058 FORMAT('0Spin orbit potential K(R)**2 FOR DISTORTED WAVE',I2)
 9503 FORMAT(18H0FORM FACTOR DATA )
 9504 FORMAT(18H PARTICLE DATA    )
 9506 FORMAT(1H0,17X,9H   LMAX =,I4,8X,6HNSTEP=,I4,8X,6HNLTR =,I4)
 9507 FORMAT(42H0QUANTUM NUMBERS FOR PROJECTILE             ,10A4)
 9508 FORMAT(1H0,17X,9H   FN   =,I4,8X,6HL    =,I4,8X,6H2*J  =,I4,8X
     1,6H2*S  =,I4,8X,6HDQ   =,F9.4)
 9906 FORMAT('0SPIN-STATISTICS DO NOT CHECK      ',6H  IS1=,I3,6H  IS2=
     1,I3,6H ISTR=,I3)
 9907 FORMAT(35H0SPIN TRANSFER DOES NOT CHECK, IS3=,I3,6H  IS4=,I3)
 9908 FORMAT(8H0WARNING,8(1H*),2X,20H0L TRANSFER EXCEEDS ,i3 )
 9999 FORMAT(I1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

      subroutine ffprint(kt, incr, ub, drf)
C
C     FORM FACTOR PRINT OUT
C
      IMPLICIT REAL*8(A-H,O-Z)
      dimension ub(800)
c
      R = DRF*float(incr)
c
      write(6,9101)
c
      MK = 2*incr - 1
      DO 3006 M=1,KT,5*incr
      MK4 = MIN0(MK+8*incr,2*KT)
      WRITE(6,9052)R ,(UB(N),UB(N+1),N=MK,MK4,2*incr) 
      MK = MK4 + 2*incr
      R=R+5.0*float(incr)*DRF
 3006 CONTINUE
      return
 9052 FORMAT(1H ,F6.2,1P10E12.4)
 9100 FORMAT(12H0FORM FACTOR,6X,9H    LTR =,I4,5X,9H   2*STR=,I4,5X
     1,                         9H   2*JTR=,I4)
 9101 FORMAT('0   R    RL,R        IM,R        RL,R+DR* N  IM,R+DR* N'
     1,             '  RL,R+DR*2N  IM,R+DR*2N  RL,R+DR*3N  IM,R+DR*3N'  
     2,             '  RL,R+DR*4N  IM,R+DR*4N' )
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      end

      SUBROUTINE FORMF(U,V,N,LAM,X, D)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(4),ZA(4),FM(4),FMA(4),RC(4),AC(4),PNLOC(4),FS(4),ECM(4)
     2,FK(4),FK2(4),ETA(4),DR(4),FMU(4),ACC(8),DELQ,AMPL(8)
     3,K,KZ,LPLUS,LPL2,IS(4),NS(4),LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,LTRP(8),JTRP(8),ISTRP(8)
     5,MLTR,IBUFF,IWORD,KMAX,NMAX,LMAX,IFLL,ICO,is_off,iextra
c
      COMMON/POTTER/DRX,AFACT(2),VFACT,SFACT,E,RM,G(4),ETAX,FKX,ETAKX
     1             ,RCX,HBARC2,ABETA(3),FLDF(3)
     2             ,NX,LAMX,KMXX,KX,IBX,LPLUSX,ICON4,NSPC,IDIRAC,ICHK
      DIMENSION U(800),V(800),X(4),D(10)
      EQUIVALENCE (G(1),FN),(G(2),FL),(G(3),FJ2),(G(4),FSS)
      data imass/ 0/
c
      ETA6=60.
      ICON4=ICON(4)
      NSPC=N
      if(imass.eq.0) then
        imass = -1
        amass = 0.0
      endif
C
C     READ IN CARD SET 5,6,OR 7   POTENTIAL CARDS
C
      READ (5,9000)E,FM(N),Z(N),FMA(N),ZA(N),RY,AC(N),PNLOC(N),FS(N),QCD
      E=E+QCD
      IF(N.EQ.2) THEN
        Q=E
        E=(ECM(1)+Q)*(FM(2)+FMA(2))/FMA(2)
      ENDIF
      IS(N)=FS(N)
      NS(N)=IS(N)+1
      DR(N)=DRF
      KMXX=KMAX
      DRX=DR(N)
      AFACT(1)=FMA(N)**.333333333
      AFACT(2)=FM (N)**.333333333
      RC(N)=ABS(RY)*AFACT(1)
      IF(RY.LT.0.0) RC(N)=RC(N)+ABS(RY)*AFACT(2)
C     CLEAR POTENTIAL STORAGE
      DO 12 M=1,KMAX
      MK=M+M-1
      U(MK  )=0.0
      U(MK+1)=0.0
      V(MK  )=0.0
      V(MK+1)=0.0
   12 CONTINUE
      RM=0.0
      IF(E.EQ.0.0) GO TO 66
C
C        ICON(10).NE.0   GIVES RELATIVISTIC KINEMATICS
C
      IF(ICON(10).NE.0) THEN
        IF(N.NE.2) GO TO 26
        IF(QCD.GT.0.0) GO TO 26
        E=E+(ECM(1)+Q)**2/(2.0*FMA(2)*AMU)
   26   CONTINUE
        elab = e
        FM1=FM(N)*AMU
        FM2=FMA(N)*AMU
        FMT=FM1+FM2
C     KLUDGE FAKE KE-LAB FOR BOUND STATES
        IF(N.GE.3) E=((E+FMT)**2-FMT**2)/(2.0*FM2)
        WLAB=E+FMT
        WCM =SQRT(2.0*E*FM2+FMT**2)
        GAMM=WLAB/WCM
        W1=GAMM*(FMT*FM1+FM2*E)/WLAB
        W2=GAMM* FM2
        ECM(N)=WCM-FMT
        IF(ICON(10).EQ.1) THEN
C       OLD KINEMATICS  SELECTED BY ICON(10) = 1
          FMU(N)=W1*W2/(W1+W2)
          VFACT=2.0*FMU(N)/HBARC**2
          SFACT=VFACT
          EFACT=VFACT
	ELSE
C       NEW KINEMATICS  SELECTED BY ICON(10) > 1
          FMU(N)=W1
          VFACT= 2.0*W1 /HBARC**2
          SFACT= 2.0*FM1/HBARC**2
          EFACT=VFACT
        ENDIF
        FMU(N)=FMU(N)/AMU
        FK2(N)=(W1**2-FM1**2)/HBARC**2
        plab = sqrt(abs((elab+fm1)**2-fm1**2))
      ELSE
C     NON RELATIVISTIC KINEMATICS
        FMU(N)=FM(N)*FMA(N)/(FM(N)+FMA(N))
        elab = e
        IF(N.LT.3) ECM(N)=E*FMU(N)/FM(N)
        IF(N.GE.3) ECM(N)=E
        VFACT=2.0*FMU(N)*AMU/HBARC**2
	SFACT=VFACT
        EFACT=VFACT
        FK2(N)=SFACT*ECM(N)
        plab = sqrt(abs(elab*2.0*amu*fm(n)))
      ENDIF
      FK(N)=SQRT(ABS(FK2(N)))
      ETAK=CHSQ*Z(N)*ZA(N)*EFACT
      ETA(N)=ETAK*0.5/FK(N)
      HBARC2=HBARC**2
C
C     ADD COULOMB AND KINETIC ENERGIES TO U
C
      RCX=RC(N)
      IF(RCX.EQ.0.0) RCX=DR(N)
      R=0.0
      FCOU=0.5*ETAK/RCX
      DO 42 M=1,KMXX
      MK=M+M-1
      R=R+DR(N)
      IF(R.GT.RCX) GO TO 40
      FC=FCOU*(3.0-(R/RCX)**2)
      GO TO 41
   40 CONTINUE
      FC=ETAK/R
   41 CONTINUE
      IF(N.LT.3) U(MK  )=U(MK  )+FK2(N)-FC
      IF(N.GE.3) U(MK+1)=U(MK+1)+FK2(N)-FC
   42 CONTINUE
      GO TO 67
   66 CONTINUE
      FK(N)=0.0
      ETA(N)=0.0
      ECM(N)=0.0
      FK2(N)=0.0
      FMU(N)=FM(N)*FMA(N)/(FM(N)+FMA(N))
      VFACT=2.0*FMU(N)*AMU/HBARC**2
      SFACT=VFACT
      ETAK=CHSQ*Z(N)*ZA(N)*VFACT
   67 CONTINUE
      IF(ICON(4).EQ.2.AND.N.GE.3) GO TO 69
      IF(N.GT.2) THEN
        Q=ECM(N)
      ELSE
        Q=ECM(N)-ECM(1)
      ENDIF
      pcm = fk(n)*hbarc
      WRITE(6,9010)N
      WRITE(6,9503)E,RY,AC(N),FS(N),plab
      WRITE(6,9504)FM(N),FMA(N),Q
      WRITE(6,9505)Z(N),ZA(N),PNLOC(N)
      WRITE(6,9500)
      RHO=FK(N)*RC(N)
      WRITE(6,9506)ECM(N),RC(N),RHO
      WRITE(6,9507)FK(N),ETA(N),DR(N),pcm
      WRITE(6,9008)
   69 CONTINUE
      FS(N)=FS(N)/2.0
      IBX=IBF(4)
      ETAX=ETA(N)
      FKX=FK(N)
      ETAKX=ETAK
      RCX=RC(N)
      LAMX=LAM
      NX=N
      LPLUSX=LPLUS
      KX=K
      ICHK=0
      IDIRAC=0
C
      CALL POTS(U,V)
C
      DR(N)=DRX
      LPLUS=LPLUSX
      IBF(4)=IBX
      IF(N.LE.2.AND.IDIRAC.NE.0.AND.PNLOC(N).EQ.0.0) PNLOC(N)=-1.0
      K=MIN0(MAX0(K,KX),KMAX)
      KZ=K
      IF(N.LE.2.OR.E.EQ.0.0) GO TO 3000
C
C     SINGLE PARTICLE ORBITAL
C
C
C     READ IN QUANTUM NUMBERS FOR SINGLE PARTICLE ORBITAL
C
      READ(5,9000)G,VTRIAL,FISW,DAMP
      ISW=FISW
      WRITE(6,9500)
      FJ0=FJ2/2.0
      FS0=FSS/2.0
      IF(VTRIAL.EQ.0.0) VTRIAL=ETA6
      WRITE(6,9508)G,VTRIAL,FISW
      IF(DAMP.EQ.0.0) GO TO 2021
      WRITE(6,9511)DAMP
 2021 CONTINUE
      FACT=(FJ0**2+FJ0-FL**2-FL-FS0**2-FS0)*0.5
c
      DO 2028 M=1,KMXX
      MK=M+M-1
      V(MK  )=U(MK  )+V(MK  )*FACT
      V(MK+1)=U(MK+1)+V(MK+1)*FACT
 2028 CONTINUE
      IF(ICON(4).NE.2) WRITE(6,9500)
      IERR=0
      IF(ICON(4).EQ.2) IERR=-1
c
      CALL BIND(V,U,DR(N),RM,FN,FL,KZ,FK(N),ETA(N),VTRIAL,ECM(N)
     1,FK2(N),ISW,IERR,D)
c
      IBF(2)=RM/DR(N)
      IBF(3)=IBF(3)+IERR
      DO 2050 M=1,KMXX
      MK=M+M-1
      V(M)=VTRIAL*V(MK  )+V(MK+1)
 2050 CONTINUE
      Anorm=1.0
      PFACT=PNLOC(N)**2/8.0
      IF(PFACT.ne.0.0) then
C
C     NON-LOCAL CORRECTION FOR SINGLE PARTICLE FUNCTION
C
      SUM=0.0
      R=0.0
      DO 2075 M=1,KZ
      MK=M+M-1
      R=R+DR(N)
      U(M)=U(M)*EXP(PFACT*(FK2(N)-V(M)))
      SUM=SUM+(U(M)*R)**2
 2075 CONTINUE
C
C     DO NOT RENORMALIZE FOR POSITIVE ENERGY
C
      IF(FK2(N).LT.0.0) then
        ANORM=1.0/SQRT(SUM*DR(N))
      ELSE
        ANORM=1.0
      ENDIF
      ENDIF
      IF(DAMP.EQ.0.0) GO TO 2080
C
C     APPLY DAMPING FACTOR, EXP(-DAMP*R) TO FORM FACTOR
C
      R=0.0
      SUM=0.0
      F1=1.0
      F2=EXP(-DAMP*DR(N))
      DO 2077 M=1,KZ
      F1=F1*F2
      R=R+DR(N)
      U(M)=U(M)*F1
      SUM=SUM+(U(M)*R)**2
 2077 CONTINUE
      Anorm=1.0/SQRT(SUM*DR(N))
 2080 CONTINUE
      DO 2100 M=KMXX,1,-1
      MK=M+M-1
      IF(M.GT.KZ) THEN
        V(M   )=0.0
        U(MK  )=0.0
      ELSE
        U(MK  )=U(M)*Anorm
      ENDIF
        U(MK+1)=0.0
 2100 CONTINUE
 3000 CONTINUE
C     Store quantum numbers for transfer back
      DO 3020 M=1,4
      X(M)=G(M)
 3020 CONTINUE
c
      RETURN
 9000 FORMAT(10F8.4)
 9008 FORMAT(21H0POTENTIAL PARAMETERS )
 9010 FORMAT( 9H0PARTICLE,I2,115(1H*))
 9500 FORMAT(1H ,3A6,5(3X,A6,F9.4))
 9503 FORMAT(18H INPUT DATA       ,9H   ELAB =,F9.4,9H   RC0  =,F9.4
     1                            ,9H   AC   =,F9.4,9H   2*STR=,F9.4
     2                            ,9h   Plab =,f9.4,9h  MeV/c  )
 9504 FORMAT(18X,9H   MASSP=,F9.4,9H   MASST=,F9.4,9H   Q    =,F9.4)
 9505 FORMAT(18X,9H   ZP   =,F9.4,9H   ZT   =,F9.4,9H   PNLOC=,F9.4)
 9506 FORMAT(18H DERIVED DATA     ,9H   ECM  =,F9.4,9H   RC   =,F9.4
     1                            ,9H   RHO  =,F9.4)
 9507 FORMAT(18X,9H   K    =,F9.4,9H   ETA  =,F9.4,9H   DR   =,F9.4
     1      ,18x,9h   P_cm =,f9.4,9h  MeV/c  )
 9508 FORMAT(18X,9H   NODES=,F9.4,9H   L    =,F9.4,9H   2*J  =,F9.4
     1          ,9H   2*S  =,F9.4,9H   VTR  =,F9.4,9H   FISW =,F9.4)
 9511 FORMAT(18X,9H   DAMP =,F9.4)
      END

