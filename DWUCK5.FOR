c$debug
      PROGRAM DWUCK5
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
      Common/n1n2ez/n1(400),n2(400),ezer(100)
      COMMON/TNT/PQ(4,8),FMTEMP(16)
C
C     STORE STANDARD ANGLE DATA
C
      ANGLE(1)=37.0
      ANGLE(2)=0.0
      ANGLE(3)=5.0
      HBARC=197.327
      AMU  =931.495
      FSCON=137.036
      CHSQ=HBARC/FSCON
c
c     kmax  = Maximum integration points for distorted waves
c     lmax  = Maximum number of partial waves
c     nmax  = Maximum number of fourier expansion coefficients
c     ifll  = Size of circular buffer in Radint
c
      KMAX =maxk
      LMAX =nls
      NMAX =imat
      IFLL =ispc1/2
      ICO  =16
      IFIRST=0
C
 1001 CONTINUE
c      CALL DW5PC  (IDAT,IFIRST)
c     CALL DW5Vax (IDAT,IFIRST)
      CALL DW5unix(IDAT,IFIRST)
c     CALL DW5ibm (IDAT,IFIRST)
      write(6,9999)ifirst
c      call ksecond(time)
c      write(*,9900)time
c      write(6,9900)time
      IFIRST=1
      ibf(3)=0
C
C     TAPE 2 STORES RADIAL INTEGRALS
C     TAPE 3 STORES LEGENDRE EXPANSION OF FORM FACTOR
C     TAPE 4 STORES DISTORTED WAVES
C     TAPE 8 STORES TNT FORM FACTOR FOR A LATER RUN OR CASE
C     TAPE 9 STORES FORM FACTORS
      REWIND 2
      REWIND 3
      REWIND 4
      REWIND 8
      REWIND 9
C
      IBF(3)=0
      WRITE(6,9502)
      WRITE(6,9804)IDAT
      icon(1)=9
C     READ CARD SET 1
      READ (5,9802,END=2)ICON,ALPHA
C     READ (5,9802)ICON,ALPHA
C
    2 continue
      WRITE(6,9803)ICON,ALPHA
      IF(ICON(1).eq.9) then
c
c      close up shop
c
        close (2)
        close (3)
        close (4)
        close (8)
        close (9)
        stop 'End of data'
      endif
C
C     THE PROGRAM MAY BE BROKEN UP INTO 3 MAIN OVERLAYS
C
C     1- ADWUCK5
C     2- BDWUCK5
C     3- CDWUCK5
C
        write(*,'('' '',15a4)')alpha
C        call ksecond(time)
C        write(*,'(a,f9.2,a)')' Adwck5 entered',time,' kseconds elapsed'
      CALL ADWCK5
      IF(IBF(3) .eq. 0) then
C        call ksecond(time)
C        write(*,'(a,f9.2,a)')' Bdwck5 entered',time,' kseconds elapsed'
        CALL BDWCK5
      ENDIF
      IF(IBF(3) .eq. 0) then
C        call ksecond(time)
C        write(*,'(a,f9.2,a)')' Cdwck5 entered',time,' kseconds elapsed'
        CALL CDWCK5
      ENDIF
c
      GO TO 1001
c
 9502 FORMAT('0CONTROL INTEGERS ')
 9802 FORMAT(20I1,15A4)
 9803 FORMAT(1H ,20I2,4X,15A4,I2)
 9804 FORMAT('  1 2 3 4 5 6 7 8 9 A B C D E F G H I J K'
     1,'    RUN IDENTIFICATION  ',44X,I4,2(1H/,I2),I4,2(1H.,I2))
 9900 FORMAT(' Elapsed time =',F9.2)
 9999 FORMAT(I1,'DWUCK5-DISTORTED WAVES U.COLORADO - PC -VAX-VERSION '
     1,' 18/Apr /1993')
      end
