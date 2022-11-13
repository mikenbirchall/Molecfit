

bool debug=1;

cpl_array* mf_derive_continuum_pars (const mf_parameters  *params, const cpl_array *fitpar, cpl_table *rangespec) {

    /* Untangle the params/fitpar arguments to access the actual continuum parameters*/
    int nmolec =     params->config->internal.molecules.n_molec;
    int nchip  =     params->config->internal.nchip;
    int ncont  = 1 + params->config->fitting.fit_continuum.n;
    int nwlc   = 1 + params->config->fitting.fit_wavelenght.n;
    int n0     = nmolec + (nwlc * nchip) + (ncont * (range - 1)); /* -> Index of where the continuum
                                                                        pars will be in the fitpar array*/

    const double *par = cpl_array_get_data_double_const(fitpar);

    /* The continuum paramters are now accessible in the par array in the subrange par[n0:n0+ncont] */

    /* Sanity Check */
    if (ncont - 1 < 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): cont_n of mf_parameters *params < 0");
    }

    /* Get the number of data points (number of lambda values)*/
    cpl_size nlam = cpl_table_get_nrow(rangespec);

    /* We will define the shift zero point of wavelength scale (wmean)  as the center of  */
    /* wavelength range. That is contiuum is defined as polynomials of (lambda-wmean)     */
    double limlam2[2] = { cpl_table_get( rangespec, MF_COL_IN_LAMBDA,  0,        NULL),
                          cpl_table_get( rangespec, MF_COL_IN_LAMBDA,  nlam - 1, NULL) };
    double wmean = (limlam2[0] + limlam2[1]) / 2;

    if(debug) cpl_msg_info(cpl_func,"MNB lam0=%f lam1=%f wmean=%f",limlam[0],limlam[1],wmean);


    /* Declare the arrays we are going to use in the form:  Ax=B */
    cpl_matrix * A   = cpl_matrix_new(nlam,ncont);
    cpl_matrix * RHS = cpl_matrix_new(nlam,1);
    cpl_matrix * SOL = cpl_matrix_new(ncont,1);
    cpl_matrix * W   = cpl_matrix_new(nlam,nlam);

    /* Initialise the Weight Matrix to 0.0. Diagonal values are defined below*/
    for(cpl_size k=0;k<nlam;k++) {
        for (cpl_size l=0;l<nlam;l++) cpl_matrix_set(W,k,l,0.0);
    }

    /* Define the RHS vector and the diagonals for the weight matrix */
    for (cpl_size k=0;k<nlam;k++)  {

        double val = cpl_table_get(rangespec, MF_COL_MOD_FLUX,  k, NULL);
        double wt1 = cpl_table_get(rangespec, MF_COL_WEIGHT,    k, NULL);
        double wt2 = cpl_table_get(rangespec, MF_COL_MOD_WEIGHT,k, NULL);
        double wt  = wt1*wt2;
        cpl_matrix_set(RHS,k,0,val);
        cpl_matrix_set(W  ,k,k,wt);

    }

    /* Initialse the SOL vector to 0.0. Uneccesary but makes me feel better*/
    for (cpl_size l=0;l<ncont;l++) cpl_matrix_set(SOL,l,0,0.0);

    /* Define the A matrix */
    for (cpl_size k=0;k<nlam;k++)  {
        double val=cpl_table_get(rangespec, MF_COL_IN_FLUX,   k, NULL);
        double lam=cpl_table_get(rangespec, MF_COL_MOD_LAMBDA,k, NULL);
        lam=lam-wmean;
        cpl_matrix_set(A,k,0,val);
        for (cpl_size l=1;l<ncont;l++) {
            val=val*lam;
            cpl_matrix_set(A,k,l,val);
        }
    }

    /* Define WA and WRHS as weighted scales of A and RHS */
    cpl_matrix* WA   =cpl_matrix_product_create(W,A);
    cpl_matrix* WRHS =cpl_matrix_product_create(W,RHS);

    /* Use SVD to solve the weighted overdetermined system */
    SOL=cpl_matrix_solve_svd(WA,WRHS);

    if(debug) cpl_matrix_dump(SOL,NULL);


    /* Make fitpar2 a copy of fitpar then update the contiuum parameters */
    cpl_array *fitpar2 = cpl_array_duplicate(fitpar);

    for  (int l1=0;l1<ncont;l1++) {
        double sval=cpl_matrix_get(SOL,l1,0);
        double fval=par[n0 + l1];
        if (debug) cpl_msg_info(cpl_func,"RANGE %lld i=%d sol=%f, cof=%f cp %f",j,l1,sval,fval,cpl_array_get_double(fitpar2,n0+l1,NULL));
        cpl_array_set_double (fitpar2,n0+l1,sval);
    }

    /* Cleanup the arrays & vectors*/
    cpl_matrix_delete(A);
    cpl_matrix_delete(W);
    cpl_matrix_delete(RHS);
    cpl_matrix_delete(WA);
    cpl_matrix_delete(WRHS);
    cpl_matrix_delete(SOL);

    return fitpar2;

}
