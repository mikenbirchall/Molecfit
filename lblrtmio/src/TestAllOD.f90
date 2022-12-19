program main

    USE Type_Kinds           , ONLY: FP, IP, DP => Double
    USE LBLRTMIO_Module

    IMPLICIT NONE

    TYPE(LBLRTM_File_type) :: ofile
    INTEGER err_stat
    INTEGER:: i, n,idx
    REAL(DP), ALLOCATABLE :: frequency(:)
    REAL(DP), ALLOCATABLE :: spectrum (:)
    REAL(DP), ALLOCATABLE :: transmission(:)
    REAL(DP) :: F1(62),F2(62),S1(62),S2(62)
    REAL(DP) :: SUM1,SUM2,PROD1,PROD2,TRANS
    INTEGER :: SIZES(62)

    CHARACTER (LEN=11) :: FNAME
    INTEGER :: NFILES

    CHARACTER(LEN=11) :: ID2ODFILENAME

    print *, 'Test lblrtmio with ODdeflt_001-ODdeflt_062 '
    print *, 'READ AND CALCULATE WITH ODdeflt files'

    FNAME=ID2ODFILENAME(1)
    CALL GET_ODFILE_DATA_SIZE(FNAME,N,F1(1),F2(1),S1(1),S2(1),.TRUE.)
    PRINT *, 'N=', N

    NFILES=62
    DO I=1,NFILES
        FNAME=ID2ODFILENAME(I)
        CALL GET_ODFILE_DATA_SIZE(FNAME,N,F1(I),F2(I),S1(I),S2(I),.FALSE.)
        SIZES(I)=N
      !  IF (I==1) CALL ParseFile(FNAME)
    ENDDO

    DO I=1,NFILES
        PRINT *, I, SIZES(I), F1(I) , '->', F2(I), ' : ', S1(I) , '->', S2(I)
    ENDDO

    PROD1=1.0
    PROD2=1.0
    DO I=1,62
        TRANS=exp(-1.0*S1(I))
        PROD1=PROD1*TRANS
        TRANS=exp(-1.0*S2(I))
        PROD2=PROD2*TRANS
    ENDDO
    PRINT *, 'PROD FOR ', F1(1), ' = ', PROD1, ' COMPARISON FROM TAPE28: 13028.00000000 ->  0.88997799 '
    PRINT *, 'PROD FOR ', F2(1), ' = ', PROD2, &
    ' COMPARISON FROM TAPE28:  3051.99984211 ->  0.73147774   AND  3052.00048828  ->  0.73157728 '

!    FNAME=ID2ODFILENAME(1)
!    CALL GET_ODFILE_DATA_SIZE(FNAME,N,.TRUE.)
!    FNAME=ID2ODFILENAME(62)
!    CALL GET_ODFILE_DATA_SIZE(FNAME,N,.TRUE.)

end program main

CHARACTER(LEN=11) FUNCTION ID2ODFILENAME(I)

    INTEGER, INTENT(IN) :: I
    CHARACTER (LEN=11) :: FNAME
    CHARACTER(LEN=1)  :: ONECHAR
    CHARACTER(LEN=2)  :: TWOCHAR

    FNAME(1:8)="ODdeflt_"
    IF (I<10) THEN
        WRITE(ONECHAR,'(I1)') I
        FNAME(9:10)="00"
        FNAME(11:11)=ONECHAR
    ELSE
        WRITE(TWOCHAR,'(I2)') I
        FNAME(9:9)="0"
        FNAME(10:11)=TWOCHAR
    ENDIF

    ID2ODFILENAME=FNAME

END FUNCTION ID2ODFILENAME

SUBROUTINE GET_ODFILE_DATA_SIZE(FNAME,N,F1,F2,S1,S2,VERBOSE)

    USE Type_Kinds           , ONLY: FP, IP, DP => Double
    USE LBLRTMIO_Module

    IMPLICIT NONE

    CHARACTER(LEN=11), INTENT(IN)  :: FNAME
    INTEGER,           INTENT(OUT) :: N
    REAL(DP), INTENT(OUT) :: F1,F2,S1,S2
    LOGICAL, INTENT(IN) :: VERBOSE

    TYPE(LBLRTM_File_type) :: ofile
    INTEGER err_stat, M
    REAL(DP), ALLOCATABLE :: frequency(:)

    err_stat = LBLRTM_File_Read(ofile, FNAME)
    IF ( err_stat /= SUCCESS ) THEN
        print *, "handle error..."
        return
    END IF
    IF (VERBOSE) THEN
        print *, "N layers=", ofile%n_Layers
        print *, "N spectra=", ofile%Layer(1)%n_Spectra
        print *, "N points=", ofile%Layer(1)%n_Points
    ENDIF

    N=ofile%Layer(1)%n_Points

    CALL LBLRTM_Layer_Frequency(oFile%Layer(1), frequency)
    m=size(frequency)
    F1=frequency(1)
    F2=frequency(m)

    S1=oFile%Layer(1)%Spectrum(1,1)
    IF (m<n ) THEN
        S2=oFile%Layer(1)%Spectrum(m,1)
        n=m
    ELSE
        S2=oFile%Layer(1)%Spectrum(n,1)
    ENDIF
    CALL LBLRTM_File_Destroy(ofile)

END SUBROUTINE GET_ODFILE_DATA_SIZE


SUBROUTINE ParseFile(FNAME)


    USE Type_Kinds           , ONLY: FP, IP, DP => Double
    USE LBLRTMIO_Module

    IMPLICIT NONE

    TYPE(LBLRTM_File_type) :: ofile
    INTEGER err_stat
    INTEGER:: i, n
    REAL(DP), ALLOCATABLE :: frequency(:)
    REAL(DP), ALLOCATABLE :: spectrum (:)
    REAL(DP), ALLOCATABLE :: transmission(:)
    REAL(DP) :: SUM1,SUM2
    CHARACTER(LEN=11) :: FNAME


    err_stat = LBLRTM_File_Read(ofile, FNAME)
    !err_stat = LBLRTM_File_Read(ofile, "TAPE3")
    !    err_stat = LBLRTM_File_Read(ofile, "TAPE12")
    IF ( err_stat /= SUCCESS ) THEN
        print *, "handle error..."
    END IF
    print *, "N layers=", ofile%n_Layers
    print *, "N spectra=", ofile%Layer(1)%n_Spectra
    print *, "N points=", ofile%Layer(1)%n_Points
    do i=1,10
        print *, i,  oFile%Layer(1)%Spectrum(i,1)
    end do

!    CALL LBLRTM_File_Inspect(ofile)

    CALL LBLRTM_Layer_Frequency(oFile%Layer(1), frequency)
    print *, size(frequency)
    print *, frequency(1), frequency(size(frequency))

    do i=1,10
        print *, i,  frequency(i), ' -> ', oFile%Layer(1)%Spectrum(i,1)
    end do
    do i=ofile%Layer(1)%n_Points-9,ofile%Layer(1)%n_Points
        print *, i,  frequency(i), ' -> ', oFile%Layer(1)%Spectrum(i,1)
    end do

    ! Now do some calculations:
    n=ofile%Layer(1)%n_Points
    allocate (spectrum(n))
    allocate (transmission(n))
    sum1=0.0d0
    sum2=0.0d0
    do i=1,ofile%Layer(1)%n_Points
        spectrum(i)=oFile%Layer(1)%Spectrum(i,1)
        sum1=sum1+frequency(i)
        sum2=sum2+spectrum(i)
        transmission(i)=exp(-1.0*spectrum(i))
    end do
    do i=n-9,n
        print *, i,  frequency(i), ' -> T -> ' , transmission(i)
    end do

    print * , 'TOTAL ', sum1,sum2

end subroutine ParseFile


