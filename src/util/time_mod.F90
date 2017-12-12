module time_mod

  implicit none

  type clock_type
    integer :: year
    integer :: month
    integer :: day
    integer :: hour
    integer :: minute
    integer :: second
    integer :: msec
    integer :: julian
    integer :: tz
    integer :: startTick
    integer :: stopTick
    integer :: currTick
    real    :: timeStep
  end type clock_type

  integer, parameter :: clock_ts_len = 28

  private

  public :: clock_type
  public :: clock_ts_len
  public :: Clock_Create, Clock_Destroy, Clock_Set, Clock_Get, Clock_Advance
  public :: IsClockDone
  public assignment(=)

  interface assignment(=)
    module procedure Clock_SetFromClock
  end interface

contains

  subroutine Clock_Create(clock)

    type(clock_type), pointer :: clock

    if (.not.associated(clock)) allocate(clock)

    clock % year   = 0
    clock % month  = 0
    clock % day    = 0
    clock % hour   = 0
    clock % minute = 0
    clock % second = 0
    clock % msec   = 0
    clock % julian = 0
    clock % tz     = 0
    clock % startTick = 0
    clock % stopTick  = 0
    clock % currTick  = 0
    clock % timeStep  = 0.
  
  end subroutine Clock_Create

  subroutine Clock_Destroy(clock)
  
    type(clock_type), pointer :: clock

    if (associated(clock)) then
      deallocate(clock)
      nullify(clock)
    end if

  end subroutine Clock_Destroy

  subroutine Clock_Set(clock, year, month, day, hour, minute, second, tz, &
                       timeStep, startTick, stopTick, currTick)

    type(clock_type), intent(inout) :: clock
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    integer, intent(in) :: minute
    integer, intent(in) :: second
    integer, optional, intent(in) :: tz
    real,    optional, intent(in) :: timeStep
    integer, optional, intent(in) :: startTick
    integer, optional, intent(in) :: stopTick
    integer, optional, intent(in) :: currTick
 
    clock % year   = year
    clock % month  = month
    clock % day    = day
    clock % hour   = hour
    clock % minute = minute
    clock % second = second
    clock % msec   = 0
    if (present(tz)) clock % tz = tz
    if (present(startTick)) then
      clock % startTick = startTick
      clock % currTick  = startTick
    end if
    if (present(stopTick))  clock % stopTick  = stopTick
    if (present(currTick))  clock % currTick  = currTick
    if (present(timeStep))  clock % timeStep  = timeStep
    
  end subroutine Clock_Set

  subroutine Clock_SetFromClock(clockOut, clockIn)
    type(clock_type), intent(inout) :: clockOut
    type(clock_type), intent(in)    :: clockIn

    clockOut % year   = clockIn % year
    clockOut % month  = clockIn % month
    clockOut % day    = clockIn % day
    clockOut % tz     = clockIn % tz
    clockOut % hour   = clockIn % hour
    clockOut % minute = clockIn % minute
    clockOut % second = clockIn % second
    clockOut % msec   = clockIn % msec
    clockOut % julian = clockIn % julian
    clockOut % startTick = clockIn % startTick
    clockOut % stopTick  = clockIn % stopTick
    clockOut % currTick  = clockIn % currTick
    clockOut % timeStep  = clockIn % timeStep
    
  end subroutine Clock_SetFromClock

  subroutine Clock_Get(clock, year, month, day, tz, hour, minute, second, msec, &
                       jday, jdate, tstamp, startTick, stopTick, currTick, timeStep)

    type(clock_type), intent(in) :: clock
    integer,          optional, intent(out) :: year
    integer,          optional, intent(out) :: month
    integer,          optional, intent(out) :: day
    integer,          optional, intent(out) :: tz
    integer,          optional, intent(out) :: hour
    integer,          optional, intent(out) :: minute
    integer,          optional, intent(out) :: second
    integer,          optional, intent(out) :: msec
    integer,          optional, intent(out) :: jday
    character(len=9), optional, intent(out) :: jdate
    character(len=clock_ts_len), optional, intent(out) :: tstamp
    integer,          optional, intent(out) :: startTick
    integer,          optional, intent(out) :: stopTick
    integer,          optional, intent(out) :: currTick
    real,             optional, intent(out) :: timeStep

    ! -- local variables
    integer :: julday

    if (present(year))   year   = clock % year
    if (present(month))  month  = clock % month
    if (present(day))    day    = clock % day
    if (present(tz))     tz     = clock % tz
    if (present(hour))   hour   = clock % hour
    if (present(minute)) minute = clock % minute
    if (present(second)) second = clock % second
    if (present(msec))   msec   = clock % msec
    if (present(jday) .or. present(jdate)) then
      julday = IW3JDN(clock % year,clock % month,clock % day) &
             - IW3JDN(clock % year,1, 1) + 1
      if (present(jday)) jday = julday
      if (present(jdate)) then
        jdate = ""
        WRITE(UNIT=jdate(1:2), FMT='(I2.2)') MOD (clock % year, 100)
        WRITE(UNIT=jdate(3:5), FMT='(I3.3)') julday
        WRITE(UNIT=jdate(6:7), FMT='(I2.2)') clock % hour
        WRITE(UNIT=jdate(8:9), FMT='(I2.2)') clock % minute
      end if
    end if
    if (present(startTick)) startTick = clock % startTick
    if (present(stopTick))  stopTick  = clock % stopTick
    if (present(currTick))  currTick  = clock % currTick
    if (present(timeStep) ) timeStep  = clock % timeStep
    if (present(tstamp)) then
      tstamp = ""
      write(tstamp, '(i4.4,2("-",i2.2),"_",2(i2.2,":"),i2.2,".",i3.3,"TZ",i2.2)') &
        clock % year, clock % month, clock % day, &
        clock % hour, clock % minute, clock % second, clock % msec, &
        clock % tz
    end if
    
  end subroutine Clock_Get

  subroutine Clock_Advance(clock, dd, hh, mm, ss, ms)

    type(clock_type), intent(inout) :: clock
    real,   optional, intent(in)    :: dd, hh, mm, ss, ms

    ! -- local variables
    integer, dimension(8) :: idat, mdat
    real,    dimension(5) :: rinc

    ! -- begin
    rinc = 0.
    if (present(dd)) rinc(1) = dd
    if (present(hh)) rinc(2) = hh
    if (present(mm)) rinc(3) = mm
    if (present(ss)) rinc(4) = ss
    if (present(ms)) rinc(5) = ms

    if (all(rinc == 0.)) then
      if (IsClockDone(clock)) return
      clock % currTick = clock % currTick + 1
      rinc(4) = clock % timeStep
    end if

    idat = (/ clock % year, clock % month, clock % day, clock % tz, &
              clock % hour, clock % minute, clock % second, clock % msec /)

    mdat = 0
    call w3movdat(rinc, idat, mdat)

    clock % year   = mdat(1)
    clock % month  = mdat(2)
    clock % day    = mdat(3)
    clock % tz     = mdat(4)
    clock % hour   = mdat(5)
    clock % minute = mdat(6)
    clock % second = mdat(7)
    clock % msec   = mdat(8)

  end subroutine Clock_Advance

  logical function IsClockDone(clock)
    type(clock_type), intent(in) :: clock

    IsClockDone = (clock % currTick > clock % stopTick)

  end function IsClockDone

#if 0
  subroutine GetJdate(yyyymmddhhmm,jdate)

    CHARACTER(len=12), INTENT(IN)  :: yyyymmddhhmm
    CHARACTER(len=9),  INTENT(OUT) :: jdate

    INTEGER year, month, day, hour, minute, jday, IW3JDN

    ! get date info from the date string
    READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
    READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
    READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
    READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour
    READ(UNIT=yyyymmddhhmm(11:12), FMT='(I2)') minute

    ! create the jdate string
    jday = IW3JDN(year,month,day) - IW3JDN(year,1, 1) + 1
    jdate = ""
    WRITE(UNIT=jdate(1:2), FMT='(I2.2)') MOD (year, 100)
    WRITE(UNIT=jdate(3:5), FMT='(I3.3)') jday
    WRITE(UNIT=jdate(6:7), FMT='(I2.2)') hour
    WRITE(UNIT=jdate(8:9), FMT='(I2.2)') minute

  end subroutine GetJdate
#endif

  INTEGER FUNCTION IW3JDN(IYEAR,MONTH,IDAY)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: IW3JDN         COMPUTE JULIAN DAY NUMBER
!   AUTHOR: JONES,R.E.       ORG: W342       DATE: 87-03-29
!
! ABSTRACT: COMPUTES JULIAN DAY NUMBER FROM YEAR (4 DIGITS), MONTH,
!   AND DAY. IW3JDN IS VALID FOR YEARS 1583 A.D. TO 3300 A.D.
!   JULIAN DAY NUMBER CAN BE USED TO COMPUTE DAY OF WEEK, DAY OF
!   YEAR, RECORD NUMBERS IN AN ARCHIVE, REPLACE DAY OF CENTURY,
!   FIND THE NUMBER OF DAYS BETWEEN TWO DATES.
!
! PROGRAM HISTORY LOG:
!   87-03-29  R.E.JONES
!   89-10-25  R.E.JONES   CONVERT TO CRAY CFT77 FORTRAN
!
! USAGE:   II = IW3JDN(IYEAR,MONTH,IDAY)
!
!   INPUT VARIABLES:
!     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
!     ------ --------- -----------------------------------------------
!     IYEAR  ARG LIST  INTEGER   YEAR           ( 4 DIGITS)
!     MONTH  ARG LIST  INTEGER   MONTH OF YEAR   (1 - 12)
!     IDAY   ARG LIST  INTEGER   DAY OF MONTH    (1 - 31)
!
!   OUTPUT VARIABLES:
!     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
!     ------ --------- -----------------------------------------------
!     IW3JDN FUNTION   INTEGER   JULIAN DAY NUMBER
!                      JAN. 1,1960 IS JULIAN DAY NUMBER 2436935
!                      JAN. 1,1987 IS JULIAN DAY NUMBER 2446797
!
!   REMARKS: JULIAN PERIOD WAS DEVISED BY JOSEPH SCALIGER IN 1582.
!     JULIAN DAY NUMBER #1 STARTED ON JAN. 1,4713 B.C. THREE MAJOR
!     CHRONOLOGICAL CYCLES BEGIN ON THE SAME DAY. A 28-YEAR SOLAR
!     CYCLE, A 19-YEAR LUNER CYCLE, A 15-YEAR INDICTION CYCLE, USED
!     IN ANCIENT ROME TO REGULATE TAXES. IT WILL TAKE 7980 YEARS
!     TO COMPLETE THE PERIOD, THE PRODUCT OF 28, 19, AND 15.
!     SCALIGER NAMED THE PERIOD, DATE, AND NUMBER AFTER HIS FATHER
!     JULIUS (NOT AFTER THE JULIAN CALENDAR). THIS SEEMS TO HAVE
!     CAUSED A LOT OF CONFUSION IN TEXT BOOKS. SCALIGER NAME IS
!     SPELLED THREE DIFFERENT WAYS. JULIAN DATE AND JULIAN DAY
!     NUMBER ARE INTERCHANGED. A JULIAN DATE IS USED BY ASTRONOMERS
!     TO COMPUTE ACCURATE TIME, IT HAS A FRACTION. WHEN TRUNCATED TO
!     AN INTEGER IT IS CALLED AN JULIAN DAY NUMBER. THIS FUNCTION
!     WAS IN A LETTER TO THE EDITOR OF THE COMMUNICATIONS OF THE ACM
!     VOLUME 11 / NUMBER 10 / OCTOBER 1968. THE JULIAN DAY NUMBER
!     CAN BE CONVERTED TO A YEAR, MONTH, DAY, DAY OF WEEK, DAY OF
!     YEAR BY CALLING SUBROUTINE W3FS26.
!
! ATTRIBUTES:
!   LANGUAGE: CRAY CFT77 FORTRAN
!   MACHINE:  CRAY Y-MP8/864, CRAY Y-MP EL2/256
!
!$$$
!
    INTEGER, INTENT(IN) :: IYEAR, MONTH, IDAY

    IW3JDN  = IDAY - 32075 &
              + 1461 * (IYEAR + 4800 + (MONTH - 14) / 12) / 4 &
              + 367 * (MONTH - 2 - (MONTH -14) / 12 * 12) / 12 &
              - 3 * ((IYEAR + 4900 + (MONTH - 14) / 12) / 100) / 4

  END FUNCTION IW3JDN


!-----------------------------------------------------------------------
  subroutine w3movdat(rinc,idat,jdat)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3MOVDAT       RETURN A DATE FROM A TIME INTERVAL AND DATE
!   AUTHOR: MARK IREDELL     ORG: WP23       DATE: 98-01-05
!
! ABSTRACT: THIS SUBPROGRAM RETURNS THE DATE AND TIME THAT IS A GIVEN
!   NCEP RELATIVE TIME INTERVAL FROM AN NCEP ABSOLUTE DATE AND TIME.
!   THE OUTPUT IS IN THE NCEP ABSOLUTE DATE AND TIME DATA STRUCTURE.
!
! PROGRAM HISTORY LOG:
!   98-01-05  MARK IREDELL
!
! USAGE:  CALL W3MOVDAT(RINC,IDAT,JDAT)
!
!   INPUT VARIABLES:
!     RINC       REAL (5) NCEP RELATIVE TIME INTERVAL
!                (DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
!     IDAT       INTEGER (8) NCEP ABSOLUTE DATE AND TIME
!                (YEAR, MONTH, DAY, TIME ZONE,
!                 HOUR, MINUTE, SECOND, MILLISECOND)
!
!   OUTPUT VARIABLES:
!     JDAT       INTEGER (8) NCEP ABSOLUTE DATE AND TIME
!                (YEAR, MONTH, DAY, TIME ZONE,
!                 HOUR, MINUTE, SECOND, MILLISECOND)
!                (JDAT IS LATER THAN IDAT IF TIME INTERVAL IS POSITIVE.)
!
! SUBPROGRAMS CALLED:
!     IW3JDN         COMPUTE JULIAN DAY NUMBER     
!     W3FS26         YEAR, MONTH, DAY FROM JULIAN DAY NUMBER
!     W3REDDAT       REDUCE A TIME INTERVAL TO A CANONICAL FORM
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      real,     intent(in) :: rinc(5)
      integer,  intent(in) :: idat(8)
      integer, intent(out) :: jdat(8)

      integer :: jldayn, jdow, jdoy
      real    :: rinc1(5),rinc2(5)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  add the interval to the input time of day and put into reduced form
!  and then compute new date using julian day arithmetic.
      rinc1(1)=rinc(1)
      rinc1(2:5)=rinc(2:5)+idat(5:8)
      call w3reddat(-1,rinc1,rinc2)
      jldayn=iw3jdn(idat(1),idat(2),idat(3))+nint(rinc2(1))
      call w3fs26(jldayn,jdat(1),jdat(2),jdat(3),jdow,jdoy)
      jdat(4)=idat(4)
      jdat(5:8)=nint(rinc2(2:5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  end subroutine w3movdat

  subroutine w3reddat(it,rinc,dinc)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3REDDAT       REDUCE A TIME INTERVAL TO A CANONICAL FORM
!   AUTHOR: MARK IREDELL     ORG: WP23       DATE: 98-01-05
!
! ABSTRACT: THIS SUBPROGRAM REDUCES AN NCEP RELATIVE TIME INTERVAL
!   INTO ONE OF SEVEN CANONICAL FORMS, DEPENDING ON THE INPUT IT VALUE.
!
!   First reduced format type (IT=-1):
!        RINC(1) is an arbitrary integer.
!        RINC(2) is an integer between 00 and 23, inclusive.
!        RINC(3) is an integer between 00 and 59, inclusive.
!        RINC(4) is an integer between 00 and 59, inclusive.
!        RINC(5) is an integer between 000 and 999, inclusive.
!      If RINC(1) is negative, then the time interval is negative.
!    
!   Second reduced format type (IT=0):
!      If the time interval is not negative, then the format is:
!        RINC(1) is zero or a positive integer. 
!        RINC(2) is an integer between 00 and 23, inclusive.
!        RINC(3) is an integer between 00 and 59, inclusive.
!        RINC(4) is an integer between 00 and 59, inclusive.
!        RINC(5) is an integer between 000 and 999, inclusive.
!      Otherwise if the time interval is negative, then the format is:
!        RINC(1) is zero or a negative integer. 
!        RINC(2) is an integer between 00 and -23, inclusive.
!        RINC(3) is an integer between 00 and -59, inclusive.
!        RINC(4) is an integer between 00 and -59, inclusive.
!        RINC(5) is an integer between 000 and -999, inclusive.
!    
!   Days format type (IT=1):
!        RINC(1) is arbitrary.
!        RINC(2) is zero.
!        RINC(3) is zero.
!        RINC(4) is zero.
!        RINC(5) is zero.
!    
!   Hours format type (IT=2):
!        RINC(1) is zero.
!        RINC(2) is arbitrary.
!        RINC(3) is zero.
!        RINC(4) is zero.
!        RINC(5) is zero.
!      (This format should not express time intervals longer than 300 years.)
!    
!   Minutes format type (IT=3):
!        RINC(1) is zero.
!        RINC(2) is zero.
!        RINC(3) is arbitrary.
!        RINC(4) is zero.
!        RINC(5) is zero.
!      (This format should not express time intervals longer than five years.)
!    
!   Seconds format type (IT=4):
!        RINC(1) is zero.
!        RINC(2) is zero.
!        RINC(3) is zero.
!        RINC(4) is arbitrary.
!        RINC(5) is zero.
!      (This format should not express time intervals longer than one month.)
!    
!   Milliseconds format type (IT=5):
!        RINC(1) is zero.
!        RINC(2) is zero.
!        RINC(3) is zero.
!        RINC(4) is zero.
!        RINC(5) is arbitrary.
!     (This format should not express time intervals longer than one hour.)
!
! PROGRAM HISTORY LOG:
!   98-01-05  MARK IREDELL
!
! USAGE:  CALL W3REDDAT(IT,RINC,DINC)
!
!   INPUT VARIABLES:
!     IT         INTEGER RELATIVE TIME INTERVAL FORMAT TYPE
!                (-1 FOR FIRST REDUCED TYPE (HOURS ALWAYS POSITIVE),
!                 0 FOR SECOND REDUCED TYPE (HOURS CAN BE NEGATIVE),
!                 1 FOR DAYS ONLY, 2 FOR HOURS ONLY, 3 FOR MINUTES ONLY,
!                 4 FOR SECONDS ONLY, 5 FOR MILLISECONDS ONLY)
!     RINC       REAL (5) NCEP RELATIVE TIME INTERVAL
!                (DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
!
!   OUTPUT VARIABLES:
!     DINC       REAL (5) NCEP RELATIVE TIME INTERVAL
!                (DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
!
! SUBPROGRAMS CALLED:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$

      integer, intent(in) :: it
      real,    intent(in) :: rinc(5)
      real,   intent(out) :: dinc(5)

!  parameters for number of units in a day
!  and number of milliseconds in a unit
!  and number of next smaller units in a unit, respectively
      integer,dimension(5),parameter:: itd=(/1,24,1440,86400,86400000/), &
          itm=(/86400000,3600000,60000,1000,1/)
!     &                                 itm=itd(5)/itd
!      integer,dimension(4),parameter:: itn=itd(2:5)/itd(1:4)
      integer,dimension(4),parameter:: itn=(/24,60,60,1000/)
      integer,parameter:: np=16
      integer :: ms
      integer iinc(4),jinc(5),kinc(5)
      real :: rp
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  first reduce to the first reduced form
      iinc=floor(rinc(1:4))
!  convert all positive fractional parts to milliseconds
!  and determine canonical milliseconds
      jinc(5)=nint(dot_product(rinc(1:4)-iinc,real(itm(1:4)))+rinc(5))
      kinc(5)=modulo(jinc(5),itn(4))
!  convert remainder to seconds and determine canonical seconds
      jinc(4)=iinc(4)+(jinc(5)-kinc(5))/itn(4)
      kinc(4)=modulo(jinc(4),itn(3))
!  convert remainder to minutes and determine canonical minutes
      jinc(3)=iinc(3)+(jinc(4)-kinc(4))/itn(3)
      kinc(3)=modulo(jinc(3),itn(2))
!  convert remainder to hours and determine canonical hours
      jinc(2)=iinc(2)+(jinc(3)-kinc(3))/itn(2)
      kinc(2)=modulo(jinc(2),itn(1))
!  convert remainder to days and compute milliseconds of the day
      kinc(1)=iinc(1)+(jinc(2)-kinc(2))/itn(1)
      ms=dot_product(kinc(2:5),itm(2:5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  next reduce to either single value canonical form
!  or to one of the two reduced forms
      if(it.ge.1.and.it.le.5) then
!  ensure that exact multiples of 1./np are expressed exactly
!  (other fractions may have precision errors)
        rp=(np*ms)/itm(it)+mod(np*ms,itm(it))/real(itm(it))
        dinc=0
        dinc(it)=real(kinc(1))*itd(it)+rp/np
      else
!  the reduced form is done except the second reduced form is modified
!  for negative time intervals with fractional days
        dinc=kinc
        if(it.eq.0.and.kinc(1).lt.0.and.ms.gt.0) then
          dinc(1)=dinc(1)+1
          dinc(2:5)=mod(ms-itm(1),itm(1:4))/itm(2:5)
        endif
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine w3reddat


  SUBROUTINE W3FS26(JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3FS26         YEAR, MONTH, DAY FROM JULIAN DAY NUMBER
!   AUTHOR: JONES,R.E.       ORG: W342       DATE: 87-03-29
!
! ABSTRACT: COMPUTES YEAR (4 DIGITS), MONTH, DAY, DAY OF WEEK, DAY
!   OF YEAR FROM JULIAN DAY NUMBER. THIS SUBROUTINE WILL WORK
!   FROM 1583 A.D. TO 3300 A.D.
!
! PROGRAM HISTORY LOG:
!   87-03-29  R.E.JONES
!   89-10-25  R.E.JONES   CONVERT TO CRAY CFT77 FORTRAN
!
! USAGE:  CALL W3FS26(JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR)
!
!   INPUT VARIABLES:
!     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
!     ------ --------- -----------------------------------------------
!     JLDAYN ARG LIST  INTEGER   JULIAN DAY NUMBER
!
!   OUTPUT VARIABLES:
!     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
!     ------ --------- -----------------------------------------------
!     IYEAR  ARG LIST  INTEGER   YEAR  (4 DIGITS)
!     MONTH  ARG LIST  INTEGER   MONTH
!     IDAY   ARG LIST  INTEGER   DAY
!     IDAYWK ARG LIST  INTEGER   DAY OF WEEK (1 IS SUNDAY, 7 IS SAT)
!     IDAYYR ARG LIST  INTEGER   DAY OF YEAR (1 TO 366)
!
!   REMARKS: A JULIAN DAY NUMBER CAN BE COMPUTED BY USING ONE OF THE
!     FOLLOWING STATEMENT FUNCTIONS. A DAY OF WEEK CAN BE COMPUTED
!     FROM THE JULIAN DAY NUMBER. A DAY OF YEAR CAN BE COMPUTED FROM
!     A JULIAN DAY NUMBER AND YEAR.
!
!      IYEAR (4 DIGITS)
!
!      JDN(IYEAR,MONTH,IDAY) = IDAY - 32075
!    &            + 1461 * (IYEAR + 4800 + (MONTH - 14) / 12) / 4
!    &            + 367 * (MONTH - 2 - (MONTH -14) / 12 * 12) / 12
!    &            - 3 * ((IYEAR + 4900 + (MONTH - 14) / 12) / 100) / 4
!
!      IYR (4 DIGITS) , IDYR(1-366) DAY OF YEAR
!
!      JULIAN(IYR,IDYR) = -31739 + 1461 * (IYR + 4799) / 4
!    &                    -3 * ((IYR + 4899) / 100) / 4 + IDYR
!
!      DAY OF WEEK FROM JULIAN DAY NUMBER, 1 IS SUNDAY, 7 IS SATURDAY.
!
!      JDAYWK(JLDAYN) = MOD((JLDAYN + 1),7) + 1
!
!      DAY OF YEAR FROM JULIAN DAY NUMBER AND 4 DIGIT YEAR.
!
!      JDAYYR(JLDAYN,IYEAR) = JLDAYN -
!     &  (-31739+1461*(IYEAR+4799)/4-3*((IYEAR+4899)/100)/4)
!
!      THE FIRST FUNCTION WAS IN A LETTER TO THE EDITOR COMMUNICATIONS
!      OF THE ACM  VOLUME 11 / NUMBER 10 / OCTOBER, 1968. THE 2ND
!      FUNCTION WAS DERIVED FROM THE FIRST. THIS SUBROUTINE WAS ALSO
!      INCLUDED IN THE SAME LETTER. JULIAN DAY NUMBER 1 IS
!      JAN 1,4713 B.C. A JULIAN DAY NUMBER CAN BE USED TO REPLACE A
!      DAY OF CENTURY, THIS WILL TAKE CARE OF THE DATE PROBLEM IN
!      THE YEAR 2000, OR REDUCE PROGRAM CHANGES TO ONE LINE CHANGE
!      OF 1900 TO 2000. JULIAN DAY NUMBERS CAN BE USED FOR FINDING
!      RECORD NUMBERS IN AN ARCHIVE OR DAY OF WEEK, OR DAY OF YEAR.
!
! ATTRIBUTES:
!   LANGUAGE: CRAY CFT77 FORTRAN
!   MACHINE:  CRAY Y-MP8/864
!
!$$$
!
    INTEGER,  INTENT(IN) :: JLDAYN
    INTEGER, INTENT(OUT) :: IYEAR,MONTH,IDAY,IDAYWK,IDAYYR 

    INTEGER :: I, J, L, N

       L      = JLDAYN + 68569
       N      = 4 * L / 146097
       L      = L - (146097 * N + 3) / 4
       I      = 4000 * (L + 1) / 1461001
       L      = L - 1461 * I / 4 + 31
       J      = 80 * L / 2447
       IDAY   = L - 2447 * J / 80
       L      = J / 11
       MONTH  = J + 2 - 12 * L
       IYEAR  = 100 * (N - 49) + I + L
       IDAYWK = MOD((JLDAYN + 1),7) + 1
       IDAYYR = JLDAYN - &
        (-31739 +1461 * (IYEAR+4799) / 4 - 3 * ((IYEAR+4899)/100)/4)

  END SUBROUTINE W3FS26

end module time_mod
