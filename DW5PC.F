      SUBROUTINE DW5PC(IDAT,ifirst)
C
c     Subroutine for opening files
c
      IMPLICIT REAL*8(A-H,O-Z)
	if(ifirst.eq.0) then
	  open( 2, FORM ='UNFORMATTED', status ='new')
	  open( 3, FORM ='UNFORMATTED', status ='new')
	  open( 4, FORM ='UNFORMATTED', status ='new')
	  open( 8, FORM ='UNFORMATTED', status ='new')
	  open( 9, FORM ='UNFORMATTED', status ='new')
C
	  CALL PCFILE(5,6,'DWUCK5 AT YOUR SERVICE   ')
C
	endif

      CALL PCDATE(IDAT)
      RETURN
      END

C
      SUBROUTINE PCFILE (INPUT, IOUTPUT, TITLE)
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*25 TITLE,INFILE,OUTFILE
      integer*2 istatus
C
C     open the input and output files
C     -------------------------------
C
	numarg = nargs() - 1
	if(numarg.ge.1) then
	  call getarg(1,infile,istatus)
	else
	  write(*,'(a\)') ' Enter INPUT  file name '
	  read (*,'(a)')  infile
	endif
	if(numarg.ge.2) then
	  call getarg(2,outfile,istatus)
	else
	  write(*,'(a\)') ' Enter OUTPUT file name '
	  read (*,'(a)') outfile
	endif

      OPEN ( INPUT,   FILE= INFILE,  STATUS='OLD' )
      OPEN ( IOUTPUT, FILE= OUTFILE, STATUS='NEW' )
C
C
C     output two line title block
C     ---------------------------
C
      WRITE (IOUTPUT,101) TITLE , INFILE 
  101 FORMAT (' ', A25 ,' **** ',' INPUT  FILE NAME = ', A )
C
      WRITE (IOUTPUT,102) OUTFILE
  102 FORMAT ( 33X, 'OUTPUT FILE NAME = ', A ,/)
      RETURN
      END


      SUBROUTINE PCDATE(IDAT)
C
c     time and date routine
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IDAT(6)
	integer*2 i1,i2,i3,i4
C                         
c	get date
c
	call getdat(i1,i2,i3)
	idat(1)=i1
	idat(2)=i2
	idat(3)=i3
c
c	get time
c
	call gettim(i1,i2,i3,i4)
	idat(4)=i1
	idat(5)=i2
	idat(6)=i3

      RETURN
      END

      subroutine second(time)
c     this subroutine interfaces MS Fortran 5.0 time routine
c     to the seconds call in the programs
c     lcsec is seconds*100
c
	real*8 time
	integer*2 ihr,imin,isec,ihun, iyr,imon,iday
	integer*4 lcini,lcsec
	logical icall
        data icall,lcini/.true., 0/
c
	call getdat(iyr,imon,iday)
	call gettim(ihr,imin,isec,ihun)
	lcsec=(3600*(24*iday + ihr) + 60*imin + isec)*100 + ihun
	if(icall) then
	  lcini = lcsec
	  icall=.false.
	endif
        time=float(lcsec - lcini)/100.
      return
      end