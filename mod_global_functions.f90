MODULE global_functions
USE global_variables
IMPLICIT NONE

CONTAINS

!..............................................................................
!
! This function computes atomic weight of a given species.
!
!..............................................................................
!
! Version 1.0 (27/02/2002)
!
!..............................................................................
!
! Input parameter(s):
!
! specs    == a name of species
!
!..............................................................................
!
! Output parameter(s):
!
! aweight  == computed atomic weight
!
!..............................................................................
REAL(wp) FUNCTION rfrac(rnum)
  IMPLICIT NONE
  INTEGER :: rnum ! The reaction number
  REAL :: reactantsum

  reactantsum = s(r(rnum)%ir1)%abundance + s(r(rnum)%ir2)%abundance
  rfrac = reactantsum/bmol
END FUNCTION rfrac 

REAL(wp) FUNCTION correct(frac)
  IMPLICIT NONE
  REAL(wp) :: frac
  REAL :: numerator,denominator

  numerator = (10**(frac*CORRFAC1)) - 1
  denominator = 100.0

  correct = numerator/denominator
END FUNCTION correct


REAL(wp) FUNCTION DILUTE(rnum)
  IMPLICIT NONE
  INTEGER :: rnum ! The reaction number
  REAL(wp) :: unity

  unity = 1.0

  dilute = correct(rfrac(rnum)) + (rfrac(rnum)**CORRFAC2)*(1.0 - correct(unity))
  IF (dilute.GT.1) THEN
    PRINT *, "ERROR: Dilution Factor greater than 1"
    PRINT *, r(rnum)%r1," + ",r(rnum)%r2," -> ",r(rnum)%p1," + ",r(rnum)%p2," + ",r(rnum)%p3
    PRINT *,"reactant fraction = ",rfrac(rnum)
    PRINT *,"Dilution factor = ",dilute
    dilute = 1.0
  ELSEIF ( dilute .LT. 0 ) THEN
    PRINT *, "ERROR: Dilution Factor less than 0"
    PRINT *, r(rnum)%r1," + ",r(rnum)%r2," -> ",r(rnum)%p1," + ",r(rnum)%p2," + ",r(rnum)%p3
    PRINT *,"reactant fraction = ",rfrac(rnum)
    PRINT *,"Dilution factor = ",dilute
    dilute = 0
  ENDIF
END FUNCTION DILUTE 


REAL(wp) FUNCTION aweight(x)
IMPLICIT NONE

! Global variable(s):
CHARACTER*10 x

! Local variable(s):
INTEGER i, first, last, n
REAL(wp) w1, w2
CHARACTER*1 s1, tmp1, multi
CHARACTER*2 s2, tmp2
LOGICAL cond
DIMENSION s1(8), w1(8), s2(7), w2(7)

PRINT *, x

! Initialization of the output:
aweight = 0.0D0

! Initialization of atom's names:
s1(1) = 'H'
s1(2) = 'C'
s1(3) = 'N'
s1(4) = 'O'
s1(5) = 'S'
s1(6) = 'P'
s1(7) = 'D'
s1(8) = 'F'
s2(1) = 'He'
s2(2) = 'Fe'
s2(3) = 'Si'
s2(4) = 'Na'
s2(5) = 'Mg'
s2(6) = 'Cl'
s2(7) = 'Xe'

! and their weights:
w1(1) = 1.00790E-00
w1(2) = 1.20110E+01
w1(3) = 1.40067E+01
w1(4) = 1.59994E+01
w1(5) = 3.20660E+01
w1(6) = 3.09738E+01
w1(7) = 2.0*w1(1)
w1(8) = 1.90000E+01
w2(1) = 4.00260E-00
w2(2) = 5.58470E+01
w2(3) = 2.80855E+01
w2(4) = 2.29898E+01
w2(5) = 2.43050E+01
w2(6) = 3.54527E+01
w2(7) = 1.31290E+02

! Initialization of borders of the species name:
!CALL len_tri2(x,10,last)
x = x(1:LEN_TRIM(x))
last = LEN_TRIM(x)


! Special species:
IF (x(1:2)=='e-' .OR. x(1:2)=='G0' .OR. x(1:2)=='G-' .OR. x(1:2)=='G+') RETURN

! First consider a special case of C10-bearing species:
IF ((x(1:3).EQ.'C10').or.(x(2:4).EQ.'C10')) THEN
	aweight = 120.11D+00
	RETURN
END IF

! Other cases:
first = 1

! Adsorbed species:
IF (x(first:first).EQ.'g') THEN
	first = 2
! Ions:
ELSE IF (x(last:last).EQ.'+') THEN
	last = last - 1
ELSE IF (x(last:last).EQ.'-') THEN
	last = last - 1
END IF

! Cyclic species
IF (x(first:first).EQ.'c') THEN
  first = first + 1
ENDIF

! Divide species onto atomic constituents:
10 CONTINUE

IF (first.GT.last) GOTO 20 ! checked the species,  exit,

! If it is not a last character in the species:
IF ((last-first).GE.1) THEN
	tmp1 = x(first:first)
    tmp2 = x(first:(first+1))
    cond = .true.

! Start a search among two-letter atoms:
    DO i = 1, 7
		IF (tmp2.EQ.s2(i)) THEN
			aweight = aweight + w2(i)
            first = first + 2
            cond = .false.
        END IF
    END DO

! Start a search among one-letter atoms:
    IF (cond) THEN
		DO i = 1, 8
			IF (tmp1.EQ.s1(i)) THEN
				multi = tmp2(2:2)
! Find a possible multiplicators:
                n = 1
                IF (multi.EQ.'2') THEN
					n = 2
                    first = first + 2
				ELSE IF (multi.EQ.'3') THEN
					n = 3
                    first = first + 2
				ELSE IF (multi.EQ.'4') THEN
                    n = 4
                    first = first + 2
				ELSE IF (multi.EQ.'5') THEN
                    n = 5
                    first = first + 2
				ELSE IF (multi.EQ.'6') THEN
					n = 6
                    first = first + 2
				ELSE IF (multi.EQ.'7') THEN
                    n = 7
                    first = first + 2
				ELSE IF (multi.EQ.'8') THEN
                    n = 8
                    first = first + 2
				ELSE IF (multi.EQ.'9') THEN
                    n = 9
                    first = first + 2

! There is no a multiplicator:
	            ELSE
                    first = first + 1
		        END IF

                aweight = aweight + n * w1(i)

            END IF

        END DO

    END IF

! Checking the last character:
ELSE

    DO i = 1, 8

		IF (x(first:last).EQ.s1(i)) THEN
			aweight = aweight + w1(i)
			first = first + 1
		END IF

    END DO

END IF

! Switch to the next character:
GOTO 10
20 CONTINUE

! Exit:
RETURN
END FUNCTION aweight

REAL(wp) FUNCTION numatoms(x)
IMPLICIT NONE

! Global variable(s):
CHARACTER*10 x

! Local variable(s):
INTEGER i, first, last, n
REAL(wp) w1, w2
CHARACTER*1 s1, tmp1, multi
CHARACTER*2 s2, tmp2
LOGICAL cond
DIMENSION s1(8), w1(8), s2(7), w2(7)

! Initialization of the output:
numatoms = 0.0d0

! Initialization of atom's names:
s1(1) = 'H'
s1(2) = 'C'
s1(3) = 'N'
s1(4) = 'O'
s1(5) = 'S'
s1(6) = 'P'
s1(7) = 'D'
s1(8) = 'F'
s2(1) = 'He'
s2(2) = 'Fe'
s2(3) = 'Si'
s2(4) = 'Na'
s2(5) = 'Mg'
s2(6) = 'Cl'
s2(7) = 'Xe'

! Initialization of borders of the species name:
!CALL len_tri2(x,10,last)
x = x(1:LEN_TRIM(x))
last = LEN_TRIM(x)

! Special species:
IF (x(1:2)=='e-' .OR. x(1:2)=='G0' .OR. x(1:2)=='G-' .OR. x(1:2)=='G+') RETURN

! First consider a special case of C10-bearing species:
IF ((x(1:3).EQ.'C10').or.(x(2:4).EQ.'C10')) THEN
	numatoms = 10
	RETURN
END IF

! Other cases:
first = 1

! Adsorbed species:
IF (x(first:first).EQ.'g') THEN
	first = 2
! Ions:
ELSE IF (x(last:last).EQ.'+') THEN
	last = last - 1
ELSE IF (x(last:last).EQ.'-') THEN
	last = last - 1
END IF

! Cyclic species
IF (x(first:first).EQ.'c') THEN
  first = first + 1
ENDIF


! Divide species onto atomic constituents:
10 CONTINUE

IF (first.GT.last) GOTO 20 ! checked the species,  exit,

! If it is not a last character in the species:
IF ((last-first).GE.1) THEN
	tmp1 = x(first:first)
    tmp2 = x(first:(first+1))
    cond = .true.

! Start a search among two-letter atoms:
    DO i = 1, 7
		IF (tmp2.EQ.s2(i)) THEN
			!aweight = aweight + w2(i)
			numatoms = numatoms + 1
            first = first + 2
            cond = .false.
        END IF
    END DO

! Start a search among one-letter atoms:
    IF (cond) THEN
		DO i = 1, 8
			IF (tmp1.EQ.s1(i)) THEN
				multi = tmp2(2:2)
! Find a possible multiplicators:
                n = 1
                IF (multi.EQ.'2') THEN
					n = 2
                    first = first + 2
				ELSE IF (multi.EQ.'3') THEN
					n = 3
                    first = first + 2
				ELSE IF (multi.EQ.'4') THEN
                    n = 4
                    first = first + 2
				ELSE IF (multi.EQ.'5') THEN
                    n = 5
                    first = first + 2
				ELSE IF (multi.EQ.'6') THEN
					n = 6
                    first = first + 2
				ELSE IF (multi.EQ.'7') THEN
                    n = 7
                    first = first + 2
				ELSE IF (multi.EQ.'8') THEN
                    n = 8
                    first = first + 2
				ELSE IF (multi.EQ.'9') THEN
                    n = 9
                    first = first + 2

! There is no a multiplicator:
	            ELSE
                    first = first + 1
		        END IF

                !aweight = aweight + n * w1(i)
                numatoms = numatoms + n

            END IF

        END DO

    END IF

! Checking the last character:
ELSE

    DO i = 1, 8

		IF (x(first:last).EQ.s1(i)) THEN
			!aweight = aweight + w1(i)
			numatoms = numatoms + 1
			first = first + 1
		END IF

    END DO

END IF

! Switch to the next character:
GOTO 10
20 CONTINUE

! Exit:
RETURN
END FUNCTION numatoms

!This function returns the number of species by its name
INTEGER FUNCTION species_idx(ss)
IMPLICIT NONE
CHARACTER*10 ss
INTEGER i

species_idx = 0

DO i = 1, nspecies
	IF (ss == s(i)%name) THEN
		IF (species_idx /= 0) THEN
			PRINT*,'ERROR: double species ',ss
			STOP
		ENDIF
		species_idx = i
	ENDIF
ENDDO

IF (species_idx>0) RETURN

IF (ss(1:LEN_TRIM(ss))/='CRP' &
    .AND. ss(1:LEN_TRIM(ss))/='CRPHOT' &
    .AND. ss(1:LEN_TRIM(ss))/='PHOTON' &
    .AND. ss(1:LEN_TRIM(ss))/='FREEZE' &
    .AND. ss(1:LEN_TRIM(ss))/='DESORB' &
    .AND. ss(1:LEN_TRIM(ss))/='IONRAD' &
    .AND. ss(1:LEN_TRIM(ss))/='PHOION' &
    .AND. ss(1:LEN_TRIM(ss))/='PHOEXC' &
    .AND. LEN_TRIM(ss)>0) THEN

    PRINT*,'WARNING: Unknown species: ',ss
    species_idx = -1
ENDIF

END FUNCTION species_idx

REAL(wp) FUNCTION delta(x,y)
IMPLICIT NONE
INTEGER :: x,y

IF (x==y) THEN
    delta = 1.0d0
ELSE
    delta = 0.0d0
ENDIF

RETURN
END FUNCTION delta

logical function EOF( unit )
integer, optional :: unit
integer :: ios

if (.not. present( unit ) ) stop "no unit given"
read( unit, *, iostat=ios )
if (ios == -1 ) then
EOF = .true.
else
EOF = .false. ; backspace( unit )
endif
end function EOF

REAL(wp) FUNCTION numelement(x,element)
  IMPLICIT NONE

  ! Global variable(s):
  CHARACTER*10 x
  CHARACTER*2 element


  ! Local variable(s):
  INTEGER i_element
  INTEGER i, first, last, n
  REAL(wp) w1, w2
  CHARACTER*1 s1, tmp1, multi
  CHARACTER*2 s2, tmp2
  LOGICAL cond
  LOGICAL twoletter
  DIMENSION s1(8), w1(8), s2(7), w2(7)


  ! Initialization of the output:
  numelement = 0.0d0

  ! Initialization of atom's names:

  s1(1) = 'H'
  s1(2) = 'C'
  s1(3) = 'N'
  s1(4) = 'O'
  s1(5) = 'S'
  s1(6) = 'P'
  s1(7) = 'D'
  s1(8) = 'F'
  s2(1) = 'He'
  s2(2) = 'Fe'
  s2(3) = 'Si'
  s2(4) = 'Na'
  s2(5) = 'Mg'
  s2(6) = 'Cl'
  s2(7) = 'Xe'

  i_element = 0
  twoletter = .FALSE.

  DO i = 1,8
    IF (TRIM(s1(i)) .EQ. TRIM(element)) THEN
      i_element = i
      twoletter = .FALSE.
      GOTO 99
    ENDIF
  ENDDO

  DO i = 1,7
    IF (TRIM(s2(i)) .EQ. TRIM(element)) THEN
      i_element = i
      twoletter = .TRUE.
      GOTO 99
    ENDIF
  ENDDO

  PRINT *, "You shouldn't be here!"
  CALL EXIT()
  99 CONTINUE

  IF ( i_element .EQ. 0 ) THEN
    PRINT *, "ERROR: Could not find element"
    RETURN
  ENDIF

  ! Initialization of borders of the species name:
  !CALL len_tri2(x,10,last)
  x = x(1:LEN_TRIM(x))
  last = LEN_TRIM(x)

  ! Special species:
  IF (x(1:2)=='e-' .OR. x(1:2)=='G0' .OR. x(1:2)=='G-' .OR. x(1:2)=='G+') RETURN

  ! First consider a special case of C10-bearing species:
  IF ((x(1:3).EQ.'C10').or.(x(2:4).EQ.'C10')) THEN
    IF ( TRIM(element).EQ.'C' ) THEN
      numelement = 10
      RETURN
    ELSE
      numelement = 0
      RETURN
    ENDIF
  END IF

  ! Other cases:
  first = 1

  ! Adsorbed species:
  IF (x(first:first).EQ.'g') THEN
    first = 2
    ! Ions:
  ELSE IF (x(last:last).EQ.'+') THEN
    last = last - 1
  ELSE IF (x(last:last).EQ.'-') THEN
    last = last - 1
  END IF

  ! Cyclic species
  IF (x(first:first).EQ.'c') THEN
    first = first + 1
  ENDIF

  ! Divide species onto atomic constituents:
  10 CONTINUE

  IF (first.GT.last) GOTO 20 ! checked the species,  exit,

  ! If it is not a last character in the species:
  IF ((last-first).GE.1) THEN
    tmp1 = x(first:first)
    tmp2 = x(first:(first+1))
    cond = .true.

    IF ( twoletter .EQV. .TRUE.) THEN
      ! Start a search among two-letter atoms:
      IF (tmp2.EQ.s2(i_element)) THEN
        !aweight = aweight + w2(i)
        numelement = numelement + 1
        first = first + 2
        cond = .false.
      END IF
      PRINT *, "You shouldn't be here!"
      CALL EXIT()
    ELSE
      ! Start a search among one-letter atoms:
      IF ((tmp1.EQ.s1(i_element)) .AND. (tmp2.NE."Si")) THEN
          multi = tmp2(2:2)
          ! Find a possible multiplicators:
          n = 1
          IF (multi.EQ.'2') THEN
            n = 2
            first = first + 2
          ELSE IF (multi.EQ.'3') THEN
            n = 3
            first = first + 2
          ELSE IF (multi.EQ.'4') THEN
            n = 4
            first = first + 2
          ELSE IF (multi.EQ.'5') THEN
            n = 5
            first = first + 2
          ELSE IF (multi.EQ.'6') THEN
            n = 6
            first = first + 2
          ELSE IF (multi.EQ.'7') THEN
            n = 7
            first = first + 2
          ELSE IF (multi.EQ.'8') THEN
            n = 8
            first = first + 2
          ELSE IF (multi.EQ.'9') THEN
            n = 9
            first = first + 2
          ELSE
            first = first + 1
          END IF

          numelement = numelement + n
      ELSE
        multi = tmp2(2:2)
        IF (multi.EQ.'2') THEN
          first = first + 2
        ELSE IF (multi.EQ.'3') THEN
          first = first + 2
        ELSE IF (multi.EQ.'4') THEN
          first = first + 2
        ELSE IF (multi.EQ.'5') THEN
          first = first + 2
        ELSE IF (multi.EQ.'6') THEN
          first = first + 2
        ELSE IF (multi.EQ.'7') THEN
          first = first + 2
        ELSE IF (multi.EQ.'8') THEN
          first = first + 2
        ELSE IF (multi.EQ.'9') THEN
          first = first + 2
        ELSE
          first = first + 1
        END IF
      END IF
    ENDIF

    ! Checking the last character:
  ELSE
    IF (x(first:last).EQ.s1(i_element)) THEN
      !aweight = aweight + w1(i)
      numelement = numelement + 1
      first = first + 1
    ELSE
      GOTO 20
    END IF
  END IF

  ! Switch to the next character:
  GOTO 10
  20 CONTINUE

  ! Exit:
  RETURN
END FUNCTION numelement

END MODULE global_functions
