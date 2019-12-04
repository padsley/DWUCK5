c$debug
      SUBROUTINE CDWCK5
c
      parameter (ispc0 = 8192, ispc1 = 8192, ispc2 = 8192, ispc3 = 2048)
      parameter (imat = 64, jmt = 2*imat, nls = 200, nl2 = 2*nls
     1, maxk = 400)
c
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (istor0 = ispc0/2, istor1 = ispc1/2)
      complex*16  D(istor0), FLL(istor1), sum(istor1)
      dimension   Plm(ispc2)
      logical D_clear,i_sym(2)
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
c
      equivalence (space0(1),D(1)),   (space1(1),FLL(1))
     1          , (space1(1),sum(1)), (space2(1),Plm(1))
c
      DATA FOURPI/12.5663706/
c
C      CALL KSECOND(TIME)
C      WRITE(6,9950)TIME
      NLTR=MLTR
      qval=ECM(2)-ECM(1)
      RMAX=FLOAT(K)*DRF 
c
c     iflag = 0 for no spin orbit, same radial elements used for all spins
      IFLAG = is(1) + is(2)
C
C       Restore spin storage
C
        IS(1)=IBF(7)
        FS(1)=FLOAT(IBF(7))/2.0
        NS(1)=IBF(7)+1
        IS(2)=IBF(8)
        FS(2)=FLOAT(IBF(8))/2.0
        NS(2)=IBF(8)+1
c       is_off .ne. 0, use full spin algebra in BETAFN routine
        is_off = is(1) + is(2) + icon(3)
c
        JR=NS(1)
        JS=NS(2)
C
      i_amp=0
      fAmass = min(fma(1), fma(2))
      fbmass = min(fm (1), fm (2))
      if    (abs(fAmass - fbmass) .lt. 0.20) then
c     special case if projectile core and target core masses are the same 
      if    (ecm(3) .gt. ecm(4)) then
        i2 = 1
      else
        i2 = 2
      endif
      elseif(abs(fma(3) - fAmass) .lt. 0.20) then
        i2 = 1
      elseif(abs(fma(3) - fbmass) .lt. 0.20) then
        i2 = 2
      endif 
c
      DO 200 II=1,NLTR
      IF(i2 .eq. 1) THEN
        LTR=LTRT(II)
        JTR=JTRT(II)
        LTP=LTRP(II)
        JTP=JTRP(II)
        IS(3)=ISTRT(II)
        IS(4)=ISTRP(II)
      ELSE
c     need to interchange form factor spins
        LTR=LTRP(II)
        JTR=JTRP(II)
        LTP=LTRT(II)
        JTP=JTRT(II)
        IS(3)=ISTRP(II)
        IS(4)=ISTRT(II)
      ENDIF
      jxr=jtr
      jxp=jtp
      istr =is(3)
      NS(3)=is(3)+1
      NS(4)=is(4)+1
C
C     Calculate normalization factor 
c     for both relativistic and nonrelativistic cases
C
      FL=1.0
      FACT=     2.0*FMU(1)/(HBARC*FK(1))**2
      FACT=FACT*2.0*FMU(2)/(HBARC*FK(2))**2
      FACT=FACT*AMU**2/FOURPI
c
      FACT=FACT*FK(2)/FK(1)
      LTMIN=MAX0(IABS(LTR-LTP),IABS(jxr-jxp)/2)+1
      LTMAX=MIN0(     LTR+LTP ,    (jxr+jxp)/2)+1
c
      DO 190 LT=LTMIN,LTMAX
      LTM=LT-1
      i_amp=i_amp+1
      JTX=JTR
      if(i_amp.eq.1.or.icon(3).eq.0) then
        D_clear = .true.
      else
        D_clear = .false.
      endif
C
      CALL BETAFN(FLL,D,LTR,LTP,JTR,JTP,LTM,jxr,jxp
     1           ,istr,iflag,is_off,D_clear,i_sym)
C
      MPLUS=JTX/2+1
      IF(ICON(8).ne.0) then
      if(icon(3).eq.0.or.(ii.eq.nltr.and.lt.eq.ltmax)) then
C
C     WRITE BETA TABLES IF ICON(8).NE.0
C
      I1 = MPLUS
      IFACT = I1*NS(1)*NS(2)
      KT=0
      IS1=-IS(1)
      DO 118 I=1,JR
      IS2=-IS(2)
      DO 116 J=1,JS
      WRITE(6,9999)ALPHA,(IDAT(ID),ID=1,3)
      WRITE(6,9905)IS2,IS1
      WRITE(6,9903)
      DO 114 LL=1,LPLUS
      LK=LL-1
      IND=KT+LK*IFACT+1
      INDF=IND+I1-1
      WRITE(6,9902)LK,(D(INDEX),INDEX=IND,INDF)
  114 CONTINUE
      KT=KT+I1
      IS2=IS2+2
  116 CONTINUE
      IS1=IS1+2
  118 CONTINUE
      endif
      endif
c
      IF(icon(3).eq.0) then
        jflag=i_amp-1
        if(jflag.ne.0.and.lt.eq.ltmax.and.ii.eq.nltr) jflag=-1
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        WRITE(6,9505)LTM,JTR,jtp,IS(3),qval,FS(1),FS(2)
c
        CALL INSIG(D,Plm,JTX,FACT,Lplus,is_off,jflag,ANGLE,NS,fm
     1            ,icon(9),icon(17),icon(19),ALPHA,IDAT,SUM)
      endif
C
  190 CONTINUE
  200 CONTINUE
c
      if(icon(3).ne.0) then
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        DO 300 II=1,NLTR
        if(fma(3).gt.fma(4)) then
          jtr=jtrt(ii)
          jtp=jtrp(ii)
          ltr=ltrt(ii)
          ltp=ltrp(ii)
        else
          jtr=jtrp(ii)
          jtp=jtrt(ii)
          ltr=ltrp(ii)
          ltp=ltrt(ii)
        endif
        LTMIN=MAX0(IABS(LTR-LTP),IABS(jxr-jxp)/2)+1
        LTMAX=MIN0(     LTR+LTP ,    (jxr+jxp)/2)+1
        DO 290 LT=LTMIN,LTMAX
        LTM=LT-1
        WRITE(6,9505)LTM,JTR,jtp,IS(3),qval,FS(1),FS(2)
  290   continue
  300   continue
        jflag=1
        CALL INSIG(D,Plm,JTX,FACT,Lplus,is_off,jflag,ANGLE,NS,fm
     1            ,icon(9),icon(17),icon(19),ALPHA,IDAT,SUM)
      endif
c
      REWIND 2
      RETURN
c
 9505 FORMAT(' ANG MOM TRANSFER ',9H   LTR  =,I4,'    2*JTR=',I4
     1,'   2*JTP=',i4,'    2*STR=',I4,'    Q    =',F9.4
     2,9H    IS1 =,F5.1,9H    IS2 =,F5.1)
 9902 FORMAT(I4,1P,10E10.3/(4X,1P,10E10.3))
 9903 FORMAT('0 L2  BETA(L2,0)          BETA(L2,1)          BETA(L2,2)'
     1,  '          BETA(L2,3)          BETA(L2,4)  ')
 9905 FORMAT('0','    MS2 =',I2,2H/2,'    MS1 =',I2,2H/2)
 9950 FORMAT('0ELAPSED TIME =',F10.3,'     SUBR. CDWUCK ENTERED')
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

      SUBROUTINE BETAFN(FLL,D,LTR,LTP,JX,JY,LT,jxr,jxp
     1                  ,istr,iflag,is_flg,D_clear,i_sym)
c
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (istor = 2**16/16)
      complex*16 FLL(istor), D(istor)
      logical D_clear,i_sym(2)
c
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(4),ZA(4),FM(4),FMA(4),RC(4),AC(4),PNLOC(4),FS(4),ECM(4)
     2,FK(4),FK2(4),ETA(4),DR(4),FMU(4),ACC(8),DELQ,AMPL(8)
     3,K,KZ,LPLUS,LPL2,IS(4),NS(4),LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,LTRP(8),JTRP(8),ISTRP(8)
     5,MLTR,IBUFF,IWORD,KMAX,NMAX,LMAX,IFLL,ICO,is_off,iextra
c
c      FLL  is radial matrix buffer
c      D    is amplitude storage
c      ltr  l of target
c      ltp  l of projectile
c      jx   j of target
c      jy   j of projectile
c      lt   l-transfer
c      iflag  = 0 use same radial matrix elements for all spin states
c      is_flg = 0 suppress spin algebra
c      D_clear  clear amplitude storage flag, used for coherent sum
c
      JR=NS(1)
      JS=NS(2)
      if(is_flg .eq. 0) then
        MPLUS = LT + 1
        JZ = 0
      else
        MPLUS = JX/2 + 1
        JZ = JY
      endif
      N = LTR+LTP-LT
      IPHAS=N-2*(N/2)
      IFACT=MPLUS*JR*JS
      INCR=LPLUS*(LT+1-IPHAS)
      J2K=(1.0+PHASEF(NS(2)))/2.0
      M2K=JX-MPLUS-MPLUS+2
      LX=LTR+LTR
      LY=LTP+LTP
      LT2=LT+LT
c
      TEMP=1.0
      TEMP1=RACAH(LY,LX,jxp,jxr,LT2,istr)
      if(is_flg .ne. 0) TEMP1=TEMP1*FLOAT(LT2+1)*PHASEF((JY+IS(3))/2)
      IF(FMA(1).GT.FMA(2)) TEMP=IBF(8)+1
      IF(FMA(1).LT.FMA(2)) TEMP=IBF(7)+1
      TEMP1=TEMP1*SQRT(TEMP*FLOAT(JZ+1))*FLOAT(LX+1)*FLOAT(LY+1)
      IF(FM(1).GT.FM(2)) TEMP1=TEMP1*PHASEF(LTR+LTP-LT)
      IF(D_clear) THEN
C
C       CLEAR AMPLITUDE STORAGE
C
        IND=LPLUS*IFACT
        DO 10 M=1,IND
        D(M)=0.0
   10   CONTINUE
      ENDIF
C
C
      IS1=-IS(1)
      DO 95 I=1,JR
      IS2=-IS(2)
      DO 90 J=1,JS
      IF(i+j.eq.2.or.iflag.ne.0) then
C
C     READ RADIAL MATRIX ELEMENTS FROM TAPE 2
C
        READ (2)(FLL(INDEX),INDEX=1,INCR)
      endif
c
      LSTOR = 0
      INDEX = 1
C
C     FINAL WAVE L LOOP
C
      DO 80 LL=1,LPLUS
      lf=LL-1
      LLX=lf+lf
      TEMP2=TEMP1*SQRT(FLOAT(LLX+1))
      JLX=LLX+IS2
      IF(JLX.LT.0) INDEX=INDEX+LT+1-IPHAS
      IF(JLX.LT.0) GO TO 75
      if(i_sym(2)) then
        if(phasef(lf).gt.0.0) then
          temp4=2.0*temp2
        else
          temp4=0.0
        endif
      else
        temp4=temp2
      endif      
      LP1=IABS(lf-LT)+1+IPHAS
      LP2=     lf+LT +1-IPHAS
      IF(lf.lt.LT) INDEX=INDEX + LP1-1-IPHAS
      IF(LP1.GT.LP2) GO TO 75
C
C     INITIAL WAVE L LOOP
C
      DO 60 LP=LP1,LP2,2
      li=lp-1
      IF (LP.le.LPLUS) then
      if(i_sym(1)) then
        if(phasef(li).gt.0.0) then
          temp3=2.0*temp4
        else
          temp3=0.0
        endif
      else
        temp3=temp4
      endif      
      LPX=LP+LP-2
      JPX=LPX+IS1
      IF(is_flg.eq.0) then
C
C     ALGEBRA FOR NO SPIN OR NO SPIN ORBIT
C
      TEMP=TEMP3*SQRT(FLOAT(LPX+1))*PHASEF(lf)
      DO 40 M=1,MPLUS
      ML=M-1
      IF(ML.le.lf) then
        MX=ML+ML
        IND=LSTOR+ML+1
        FACT=VCC(LLX,LPX,LT2,MX,0)
     1   *SQRT(YXFCT(lf+ML,lf-ML))*TEMP
        D(IND  )=D(IND  )+FLL(INDEX  )*FACT
      endif
   40 CONTINUE
      else
C
C     ALGEBRA FOR SPIN WITH SPIN ORBIT or coherent sum
C
      IF(JPX.ge.0) then
      TEMP=TEMP3*SQRT(FLOAT(LPX+1)*FLOAT(JLX+1))
     1          *WINEJ(LLX,LT2,LPX,IS(2),JY,IS(1),JLX,JX,JPX)
      KT=0
      MSP=-IS(1)
      DO 57 MS1=1,JR
      MS =-IS(2)
      DO 55 MS2=1,JS
      VCP=VCC(LPX,IS(1),JPX,0,MSP)
      DO 50 M=1,MPLUS
      MK=M+M-1
      MX=MK-1+M2K
      ML2=MSP-MX-MS
      ML=IABS(ML2/2)
      IF(ML.le.lf) then
        IND=LSTOR+KT+M
        FACT=VCP*VCC(JLX,JX,JPX,MSP-MX,MX)*VCC(LLX,IS(2),JLX,ML2,MS)
     1   *SQRT(YXFCT(lf+ML,lf-ML))*TEMP
        D(IND  )=D(IND  )+FLL(INDEX  )*FACT
      endif
   50 CONTINUE
      KT=KT+MPLUS
      MS =MS +2
   55 CONTINUE
      MSP=MSP+2
   57 CONTINUE
      endif
      endif
      endif
   59 CONTINUE
      INDEX=INDEX+1
   60 CONTINUE
   75 CONTINUE
      LSTOR=LSTOR+IFACT
   80 CONTINUE
      IS2=IS2+2
   90 CONTINUE
      IS1=IS1+2
   95 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE INSIG(D,Plm,JTR,FACTR,Lplus,is_flg,JFLAG, ANGLE, NS, fm
     1         , icon9, icon17, icon19, ALPHA, IDAT, SUM )
c
c     Calculate inelastic cross sections and spin observables.
c
c     D()     Final state amplitudes
c     Plm()   Legendre polynomial storage = Lplus*(Jtr+1/2)
c     JTR     2*J_transfer
c     FACTR   normalization factor from CDWUCK
c     Lplus   number of partial waves, Lmax+1
c     is_flag
c     JFLAG
c     ANGLE() N_angles, Angle_0, Del_angle
c     NS()    spin multiplicities for the distorted waves
c     fm()    atomic masses of projectiles 
c     icon9   .ne. 0 enables inelastic observable printing to output file
c     icon17  .ne. 0 enables output of spin observables to Fortran files 20, 21
c     icon19  .ne. 0 disables all printinting to output file
c     ALPHA   Run identification
c     IDAT()  date and time,  (yy,mm,dd)  (hh:mm:ss)
c     SUM()   angle dependent amplitude storage = JR*JS*Lplus*(Jtr+1/2)
c     
c***********************************************************************
c
c     nx     maximum number of angles for output
c
      parameter (nx = 200, npol = 10)
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical i_open20, i_out20, i_20flag
c
      double complex D(*), SUM(*), SUM1
      dimension sigplt(nx), asyplt(nx), sigsum(nx), IDAT(6)
     1   , Polz(npol), plm(*), ANGLE(5), fm(2), NS(4), IS(4)
     2   ,ALPHA(15)
c
      equivalence (SIG, Polz(1))

      DATA i_20flag / .true./
C
c     --------------------------------------------------------
c       Calculate cross section and spin observables and
c       output to disk file 20 if( icon(17) .ne. 0)
c
      
      if(icon17 .eq. 2 .or. icon17 .eq. 3) then
        i_out20  = .true.
        i_open20 = .true.
        if(i_20flag) then
          i_open20 = .true.
        endif
      else    
        i_open20 = .false.
        i_out20  = .false.
c
      endif
      NTHETA = ANGLE(1)
      NTHETA = MIN0(NTHETA,NX)
      THETA1 = ANGLE(2)
      DTHETA = ANGLE(3)
c
      JR=NS(1)
      JS=NS(2)
c     IS are 2*spin of the distorted waves      
      IS(1) = NS(1) - 1
      IS(2) = NS(2) - 1
      IS(3) = NS(3) - 1
      if(is(1).eq.2.and.fm(1).eq.0.0) then
c     initial state average factor for Gamma ray initial state
        FACTA=sqrt(FACTR*3.0/2.0)
      else
        FACTA=sqrt(FACTR)
      endif
c      
      M2K=(1.0-PHASEF(IS(3)))/2.0
      IF(is_flg.EQ.0) M2K=0
      NPLUS=(JTR+IS(1)+IS(2))/2+1
      MPLUS=JTR/2+1
      IFACT = MPLUS*JR*JS
      IF(NTHETA.EQ.0) GO TO 230
c     write column label for observable output
      if(icon19 .eq. 0) then
       WRITE(6,9000)         
      endif
      TotSig=0.0
c
c     clear sigsum storage 
c
      IF(JFLAG.eq.0) then
        DO 5 NTH=1,NTHETA
        SIGSUM(NTH)=0.0
    5   CONTINUE
      endif
C
      THETA = THETA1
      DO 110 NTH=1,NTHETA
      CALL LGNDR(Plm,NPLUS,LPLUS,THETA)
c
      Index1 = JS*JR*((JTR+1)/2)
      Index2 = Index1+1
      DO 100 M=1,MPLUS
      M2 = 2*(M-1)+M2K
      KT=M
      IS1=-IS(1)
c     Loop for initial spin states
      DO 70 I=1,JR
      IS2=-IS(2)
c     Loop for final   spin states
      DO 60 J=1,JS
      ML=-(M2+IS2-IS1)/2 
c      
      if(ML .lt. 0)then
        PHAS1=PHASEF(ML)
      else
        PHAS1=1.0
      endif  
      if(ML .gt. 0)then
        PHAS2=PHASEF(ML)
      else
        PHAS2=1.0
      endif             
c      
      ML1=IABS(ML)*LPLUS
      SUM1=0.0
c
      IND=KT
      DO 40 LL=1,LPLUS
      ML1=ML1+1
      SUM1 = SUM1+D(IND)*Plm(ML1)
C
C     CALCULATE TOTAL INELASTIC SIGMA
C
      IF(NTH .eq. 1) THEN
      L=LL-1
      MM = iabs(ML) 
        if(MM .le. L) then
          FACT=conjg(D(IND))*D(IND)*YXFCT(L-MM,L+MM)/FLOAT(2*L+1)
          IF(M2.NE.0) FACT=FACT*2.0
          TotSig=TotSig+FACT
        endif
      ENDIF
      IND=IND+IFACT
   40 CONTINUE
        Index1 = Index1+1
        SUM(Index1) = SUM1*PHAS1 *FACTA
      if(M2 .ne. 0) then
        Index2 = Index2-1
        SUM(Index2) = SUM1*PHAS2 *FACTA
      endif
c	if(nth.eq.1) write(*,'(a,4i3, 1p4e12.4)')
c     1  ' Is2,Is1 M, ML :',is2,is1,M,ML,SUM(Index1),SUM(Index2)
      KT=KT+MPLUS
      IS2=IS2+2
   60 CONTINUE
      IS1=IS1+2
   70 CONTINUE
c
  100 CONTINUE
c
      Maxi = JTR + 1
      Max1 = 1
      CALL POLFCT(Max1,Maxi,JR,JS,Theta,Polz,SUM
     1  ,i_open20,i_out20,nth,ntheta,ALPHA,IDAT)
c
      SIGSUM(NTH)=SIGSUM(NTH)+SIG
      WRITE(6,9001)THETA, (Polz(I),I=1,10),theta
        SIGPLT(NTH)=SIG
        ASYPLT(NTH)=Polz(3)
      THETA=THETA+DTHETA
  110 CONTINUE
c
      TotSig =TotSig *12.5663706*FACTA**2/float(JR)
c      
      if(i_out20) then
        WRITE(20,9905)' Totl_Sig',TotSig
      endif
c
      if(icon19 .eq. 0) then
c
      WRITE(6,9002)TotSig
c      
      IF(NTHETA.ne.0.and.ICON9.ne.0) then
      WRITE(6,9011)ALPHA,(IDAT(I),I=1,3)
      NTH=NTHETA
      CALL DWPLOT(NTH,ICON9,SIGPLT,DTHETA,THETA1)
      NTEMP=-10
      IF(IS(1).NE.0) THEN
        WRITE(6,9013)ALPHA,(IDAT(I),I=1,3)
        CALL DWPLOT(NTH,NTEMP  ,ASYPLT,DTHETA,THETA1)
      ENDIF
C
      IF(JFLAG.EQ.-1) THEN
C
C     PRINT SUMMED INCOHERENT CROSS SECTIONS
C
      WRITE(6,9011)ALPHA,(IDAT(I),I=1,3)
      WRITE(6,9004)
      THETA=THETA1
      DO 220 NTH=1,NTHETA
      SIG1=SIGSUM(NTH)
      WRITE(6,9001)THETA,SIG1
      THETA=THETA+DTHETA
  220 CONTINUE
      if(ntheta.ne.0.and.icon9.ne.0) then
        NTH=NTHETA
        WRITE(6,9011)ALPHA,(IDAT(I),I=1,3)
        CALL DWPLOT(NTH,ICON9,SIGSUM,DTHETA,THETA1)
      endif
      ENDIF
      ENDIF
      endif
  230 continue
      RETURN
c
 9000 FORMAT('0 Theta',' Inelsig,fm**2'
     1,   '      Polz','       Asy','       Ayy'  
     2,   '       A22','       A21','       A20'
     3,   '       T22','       T21','       T20'
     4,   '  Theta')
 9001 FORMAT(0PF7.2,1PE14.4,0P9F10.4,0pf7.2)
 9002 FORMAT('0Tot-sig',1PE13.4)
 9004 FORMAT('0 Incoherent summed sigma')
 9011 FORMAT(12H1Inelastic  ,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
 9013 FORMAT(12H1Asymmetry  ,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
 9905 FORMAT(A9,1pe13.4)
       END
