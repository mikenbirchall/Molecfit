cpl_boolean mf_QWKAParseEnvVar(QWKA_FLAG FLAG) {

    /* Default Flag values*/
    cpl_boolean VERBOSE_FLAG   = CPL_FALSE;
    cpl_boolean COMPARE_FLAG   = CPL_FALSE;
    cpl_boolean TERMINATE_FLAG = CPL_FALSE;
    cpl_boolean TAPE5_FLAG     = CPL_FALSE;

    /* Get the env var*/
    char* FLAG_STR=NULL;
    FLAG_STR=getenv("QWKMFA");

    if (FLAG_STR!=NULL) {

        /* Parse the env var string into flag states */
        for (size_t i=0; FLAG_STR[i]!='\0'; i++) {
            if (FLAG_STR[i]=='V')   VERBOSE_FLAG = CPL_TRUE;
            if (FLAG_STR[i]=='v')   VERBOSE_FLAG = CPL_FALSE;
            if (FLAG_STR[i]=='C')   COMPARE_FLAG = CPL_TRUE;
            if (FLAG_STR[i]=='c')   COMPARE_FLAG = CPL_FALSE;
            if (FLAG_STR[i]=='T') TERMINATE_FLAG = CPL_TRUE;
            if (FLAG_STR[i]=='t') TERMINATE_FLAG = CPL_FALSE;
            if (FLAG_STR[i]=='F')     TAPE5_FLAG = CPL_TRUE;
            if (FLAG_STR[i]=='f')     TAPE5_FLAG = CPL_FALSE;
        }

}

    /* Return the requred flag value*/
    switch (FLAG) {
        case VERBOSE  : return VERBOSE_FLAG;
        case COMPARE  : return COMPARE_FLAG;
        case TERMINATE: return TERMINATE_FLAG;
        case TAPE5    : return TAPE5_FLAG;
        default       : cpl_msg_info(cpl_func,"Error bad Flag");
    }

    return CPL_FALSE;
}

