program main

    USE LBLRTMIO_Module
    IMPLICIT NONE

    TYPE(LBLRTM_File_type) :: ofile
    INTEGER err_stat
    print *, 'Test lblrtmio'

    err_stat = LBLRTM_File_Read(ofile, "ODint_062")
    IF ( err_stat /= SUCCESS ) THEN
        print *, "handle error..."
    END IF

    CALL LBLRTM_File_Inspect(ofile)


end program main
