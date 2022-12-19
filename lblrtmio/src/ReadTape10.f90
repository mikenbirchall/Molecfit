program main

    USE Type_Kinds           , ONLY: FP, IP, DP => Double
    USE LBLRTMIO_Module
    IMPLICIT NONE

    TYPE(LBLRTM_File_type) :: ofile
    INTEGER err_stat

    REAL(DP), ALLOCATABLE :: frequency(:)
    REAL(DP), ALLOCATABLE :: spectrum (:)
    REAL(DP), ALLOCATABLE :: transmission(:)

    INTEGER:: I, NPTS


    print *, 'READ TAPE10'

    err_stat = LBLRTM_File_Read(ofile, "TAPE10")
    IF ( err_stat /= SUCCESS ) THEN
        print *, "handle error..."
    END IF

    !CALL LBLRTM_File_Inspect(ofile)

    print *, "N layers=", ofile%n_Layers
    print *, "N spectra=", ofile%Layer(1)%n_Spectra
    print *, "N points=", ofile%Layer(1)%n_Points


    CALL LBLRTM_Layer_Frequency(oFile%Layer(1), frequency)


    print *, size(frequency)
    print *, frequency(1), frequency(size(frequency))


    do i=1,10
        print *, i,  frequency(i), ' -> ', oFile%Layer(1)%Spectrum(i,1)
    end do
    do i=ofile%Layer(1)%n_Points-9,ofile%Layer(1)%n_Points
        print *, i,  frequency(i), ' -> ', oFile%Layer(1)%Spectrum(i,1)
    end do

    npts=ofile%Layer(1)%n_Points
    allocate(spectrum(npts))
    allocate(transmission(npts))
    do i=1,npts
        spectrum(i)=oFile%Layer(1)%Spectrum(i,1)
        transmission(i)=exp(-1.0*spectrum(i))
    end do

    do i=1,10
        print *, i,  frequency(i), ' -> ', spectrum(i), ' -> ', transmission(i)
    end do
    do i=npts-9,npts
        print *, i,  frequency(i), ' -> ', spectrum(i), ' -> ', transmission(i)
    end do





end program main
