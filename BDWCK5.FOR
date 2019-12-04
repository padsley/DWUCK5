c$debug
      SUBROUTINE BDWCK5
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical post_strip
c
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(4),ZA(4),FM(4),FMA(4),RC(4),AC(4),PNLOC(4),FS(4),ECM(4)
     2,FK(4),FK2(4),ETA(4),DR(4),FMU(4),ACC(8),DELQ,AMPL(8)
     3,K,KZ,LPLUS,LPL2,IS(4),NS(4),LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,LTRP(8),JTRP(8),ISTRP(8)
     5,MLTR,IBUFF,IWORD,KMAX,NMAX,LMAX,IFLL,ICO,is_off,iextra
c
      COMMON/BDWCK_5/LTR,LTP,JTR,JTP,T(4),QMAX(2),DP,PN,DQ,QN
      COMMON/TNT    /PQ(4,8),FMTEMP(16)
c
      DIMENSION AKEX(2)
c
      EQUIVALENCE (ACC(7),AKEX(1))
c
c      Test for consistent transferred particle masses
      if(icon(11) .eq. 0 .and. abs(FM(3)-FM(4)) .gt. 0.20) go to 600
C
      IF(FMA(1).LT.FMA(2)) THEN
C     HERE FOR (D,P)    A(a,b)B  B = A+x, a = b+x
        I1=2
        I2=1
c		First   form factor,  k_d - (A/B)k_p
        T(1)=1.0
        T(2)=FMA(I2)/FMA(I1)
c		Second  form factor,  (a/b)k_d - k_p
        T(3)=FM(I1)/FM(I2)
        T(4)=1.0
C     NO RECOIL OPTION
      IF(ICON(17).EQ.1) T(4)=T(3)*T(2)
      ELSE
C     HERE FOR (P,D)    B(b,a)A  B = A+x, a = b+x
        I1=1
        I2=2
c		First   form factor,  (A/B)k_p - k_d
        T(1)=FMA(I2)/FMA(I1)
        T(2)=1.0
c		Second  form factor,   k_p - (a/b)k_d
        T(3)=1.0
        T(4)=FM(I1)/FM(I2)
C     NO RECOIL OPTION
      IF(ICON(17).EQ.1) T(3)=T(4)*T(1)
      ENDIF
C
c     Fill in for Target mass 4 for icon(11) .ne. 0
c
      if(icon(11) .ne. 0) then
        if(abs(FMA(3)-FMA(i2)) .lt. 0.2) then
          FMA(4) = FM(i1)
        else
          FMA(4) = FMA(i2)
        endif
      endif
C
      IF(ICON(13).NE.0) WRITE(6,9103)
c
      IF    (ABS(FMA(3)-FMA(I2)) .lt. 0.20) then
        post_strip = .true.
c       Special case when projectile core mass is same as target core mass
        if(ABS(FMA(3)-FMA(4)) .lt. 0.20) then
        if    (FMA(1) .lt. FMA(2) .and. ECM(3) .lt. ECM(4)) then
          post_strip = .false.
        elseif(FMA(1) .gt. FMA(2) .and. ECM(3) .gt. ECM(4)) then
          post_strip = .false.
        endif
        endif
      ELSEIF(ABS(FMA(4)-FMA(I2)) .lt. 0.20) then
        post_strip = .false.
        if(ABS(FMA(3)-FMA(4)) .lt. 0.20) then
        if    (FMA(1) .lt. FMA(2) .and. ECM(3) .gt. ECM(4)) then
          post_strip = .true.
        elseif(FMA(1) .gt. FMA(2) .and. ECM(3) .lt. ECM(4)) then
          post_strip = .true.
        endif
        endif
      ELSE
        go to 500
      ENDIF
      if(post_strip) then
C
C       HERE IF particle 3 IS FOR TARGET
C       HERE IF particle 4 IS FOR PROJECTILE
c
c       Post  stripping and prior pick up form
C
        IF(ICON(13).EQ.0) THEN
          QMAX(1)=AKEX(I2)+FMA(I2)*AKEX(I1)/FMA(I1)
          QMAX(2)=FM(I1)*AKEX(I2)/FM(I2)+AKEX(I1)
        ELSE
          QMAX(1)=SQRT(AKEX(I2)**2+(FMA(I2)*AKEX(I1)/FMA(I1))**2)
          QMAX(2)=SQRT((FM(I1)*AKEX(I2)/FM(I2))**2+AKEX(I1)**2)
        ENDIF
c
      else
C
C       HERE IF particle 3 IS FOR PROJECTILE
C       HERE IF particle 4 IS FOR TARGET
c
c       Prior stripping and post  pick up form
C
        IF(ICON(13).EQ.0) THEN
          QMAX(1)=FM(I1)*AKEX(I2)/FM(I2)+AKEX(I1)
          QMAX(2)=AKEX(I2)+FMA(I2)*AKEX(I1)/FMA(I1)
        ELSE
          QMAX(1)=SQRT((FM(I1)*AKEX(I2)/FM(I2))**2+AKEX(I1)**2)
          QMAX(2)=SQRT(AKEX(I2)**2+(FMA(I2)*AKEX(I1)/FMA(I1))**2)
        ENDIF
C       interchange momenta  between form factors
        TEMP=T(1)
        T(1)=T(3)
        T(3)=TEMP
        TEMP=T(2)
        T(2)=T(4)
        T(4)=TEMP
      endif
C
      NLTR=MLTR
      DO 200 II=1,NLTR
      LTR=LTRT(II)
      LTP=LTRP(II)
      JTR=JTRT(II)
      JTP=JTRP(II)
      IS(3)=ISTRT(II)
      DP=PQ(1,II)
      PN=PQ(2,II)
      DQ=PQ(3,II)
      QN=PQ(4,II)
C
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
c	write(6,'(a,4f10.5)')'  T(i) ',T
      IF(ICON(2).eq.0) then
C
C     This section for one nucleon-like transfer
C
      CALL BADWK5
C     CALL EXPAND(LTR,DQR,QMAX)
C     CALL XFORM(LTR,LTP,T,DQ)
      else
C
C     This section for two nucleon transfer
C
      CALL BXDWK5
C     CALL EXPAN2
C     CALL DSTRI2
C     CALL XFORM2(LTR,LTP,T,DQR)
C
      endif
      REWIND 3
      CALL BBDWK5
C     CALL RADINT(LTR,LTP,JTR,JTP,T)
C
      REWIND 3
  200 CONTINUE
      END FILE 2
      REWIND 2
      REWIND 9
      RETURN
  500 CONTINUE
      WRITE(6,9100)
      IBF(3) = 1
      RETURN
  600 continue
      write(6,9101)FM(3),FM(4)
      IBF(3) = 1
      return
c
 9100 FORMAT('0MASSES ARE NOT CONSISTENT. EXIT FROM Subroutine BDWUCK ')
 9101 format('0Transferred particle masses are not the same,  FM(3) ='
     1  ,f8.2,'     FM(4) =',f8.2 )
 9103 FORMAT('0MOMENTUM SPACE TRUNCATED TO ANGLES.LE.90 DEGREES ')
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END
