/*
 * This file is part of the ESO Telluric Correction Library
 * Copyright (C) 2001-2018 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include "mf_kernel_synthetic.h"
#include "mf_kernel_user.h"

#include "mf_convolution.h"

/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Enumeration types
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Defines
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Global variables
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Macros
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Structured types
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/*  */
static cpl_error_code mf_convolution_rebin(
    cpl_table                *outspec,
    const char               *outlam,
    const char               *outflux,
    const cpl_table          *inspec,
    const char               *inlam,
    const char               *influx);

/*  */
static cpl_error_code mf_convolution_mod_wave_grid(
    cpl_table                *spec,
    const mf_parameters      *params,
    const cpl_array          *fitpar,
    const int                chip);

/*  */
static cpl_error_code mf_convolution_telback(
    cpl_table                *spec,
    const mf_parameters      *params,
    const cpl_array          *fitpar);

/*  */
static cpl_error_code mf_convolution_mod_continuum(
    cpl_table                *spec,
    const mf_parameters      *params,
    const cpl_array          *fitpar,
    const int                range);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_convolution       .
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief Execute convolution
 *
 * @param params             .
 * @param correct_spectrum   .
 * @param last_call          .
 * @param spec_out           .
 * @param range_status       .
 * @param fitpar             .
 * @param mpfit_calls        .
 * @param spec               .
 * @param kernel             .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_convolution(
    mf_parameters      *params,
    const cpl_boolean        correct_spectrum,
    const cpl_boolean        last_call,
    cpl_table                **spec_out,
    cpl_error_code           *range_status,
    cpl_array                *fitpar,
    int                      mpfit_calls,
    cpl_table                *spec,
    cpl_matrix               *kernel)
{
    /*** Second Part: Convolution ***/

    /* Set model columns to default values */
    cpl_size m = cpl_table_get_nrow(spec);
    cpl_table_fill_column_window(spec, MF_COL_MOD_LAMBDA, 0, m, 0.);
    cpl_table_fill_column_window(spec, MF_COL_MOD_SCALE,  0, m, 1.);
    cpl_table_fill_column_window(spec, MF_COL_MOD_FLUX,   0, m, 0.);
    cpl_table_fill_column_window(spec, MF_COL_MOD_WEIGHT, 0, m, 0.);
    cpl_table_fill_column_window(spec, MF_COL_DEV,        0, m, 0.);

     /* MNB*/
     cpl_table* flux0 = cpl_table_new(m);
     cpl_table_duplicate_column(flux0,"FLUX0",spec,MF_COL_DEV);
     cpl_boolean MFIT_SVD = CPL_FALSE;
     char * envStr=NULL;
     envStr=getenv("MFIT_SVD");
     if (envStr!=NULL) {
         if (envStr[0]=='1') {
             cpl_msg_info(cpl_func,"MFIT_SVD==TRUE");
             MFIT_SVD=CPL_TRUE;
         } else {
            cpl_msg_info(cpl_func,"MFIT_SVD==FALSE");
         }
     }
     if (MFIT_SVD)             cpl_msg_info(cpl_func,"Using SVD");



    /* Calculate one molecular spectrum for full/selected wavelength range? */
    cpl_boolean single_spectrum = params->config->internal.single_spectrum;

    /* Adapt model spectrum for each part (chip) of the observed spectrum */
    int nrange = params->config->internal.n_range;
    cpl_table *modspec = NULL;
    for (cpl_size j = 0; j < nrange; j++) {

        /* Skip empty ranges */
        if (cpl_table_get(params->rangetab, MF_COL_WN_END, j, NULL) != 0.) {

            cpl_error_code codestat = CPL_ERROR_NONE;

            /* Convert wavenumbers to wavelengths in the RANGE out spectrum of LBLRTM */
            /*cpl_table *modspec = NULL;*/
            /*modspec  = cpl_table_new(0);*/
            if (!single_spectrum || (single_spectrum && j == 0) ) {

                modspec  = cpl_table_new(0);
                codestat = range_status[j];

                if (codestat == CPL_ERROR_NONE && !(spec_out[j]) ) {

                    codestat = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                                 "LBLRTM spectrum range=%lld empty!", j);
                } else {

                    cpl_table_set_size(        modspec, cpl_table_get_nrow(spec_out[j]));
                    cpl_table_new_column(      modspec, MF_COL_IN_LAMBDA, CPL_TYPE_DOUBLE);
                    cpl_table_new_column(      modspec, MF_COL_IN_FLUX,   CPL_TYPE_DOUBLE);
                    cpl_table_copy_data_double(modspec, MF_COL_IN_LAMBDA, cpl_table_get_data_double(spec_out[j], MF_COL_IN_LAMBDA));
                    cpl_table_copy_data_double(modspec, MF_COL_IN_FLUX,   cpl_table_get_data_double(spec_out[j], MF_COL_IN_FLUX  ));

                    /* Convert fluxes if radiance spectrum */
                    const double conv = 1e-4 / (MF_CONVOLUTION_LAM_UNIT * CPL_PHYS_C * CPL_PHYS_H * MF_CONVOLUTION_SR_IN_ARCSEC2);
                    if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE && correct_spectrum) {
                        cpl_table_divide_columns( modspec, MF_COL_IN_FLUX, MF_COL_IN_LAMBDA);
                        cpl_table_multiply_scalar(modspec, MF_COL_IN_FLUX, conv            );
                    }

                    /* Avoid "nan" (not a number) in input table */
                    cpl_boolean isnanum = CPL_FALSE;
                    for (cpl_size i = 0; i < cpl_table_get_nrow(modspec); i++) {
                        if (isnan(cpl_table_get(modspec, MF_COL_IN_FLUX, i, NULL)) != 0) {
                            cpl_table_set(modspec, MF_COL_IN_FLUX, i, -9e99);
                            isnanum = CPL_TRUE;
                        }
                    }

                    /* Print warning message in the case of "nan" values */
                    if (isnanum == CPL_TRUE) {
                        cpl_msg_warning(cpl_func, "(mf_convolutio) NaN values -> weight = 0");
                    }
                }
            };/*(!single_spectrum || (single_spectrum && j == 0) ) */
            
            if (codestat != CPL_ERROR_NONE) cpl_msg_info(cpl_func,"LBLRTM ERROR  so returning err");
            if (!modspec) cpl_msg_info(cpl_func,"!not modspec so returning err");
            if (!modspec) return CPL_ERROR_ILLEGAL_INPUT;

            /* LBLRTM errors -> output: model spectrum = observed spectrum */
            if (codestat != CPL_ERROR_NONE) {
                cpl_table_delete(modspec);
                return codestat;
            }

            /* Extract wavelength grid of observed spectrum for each range */
            cpl_table_unselect_all(spec);

            cpl_size  nsel         = cpl_table_or_selected_int(spec, MF_COL_MOD_RANGE, CPL_EQUAL_TO, j + 1);
            cpl_table *rangespec   = cpl_table_extract_selected(spec);
            cpl_array *origselrows = cpl_table_where_selected(spec);
            cpl_array *selrows     = cpl_array_cast(origselrows, CPL_TYPE_INT);

            cpl_array_delete(origselrows);
            cpl_table_select_all(spec);

            cpl_table_erase_column(rangespec, MF_COL_CHIP);
            cpl_table_name_column( rangespec, MF_COL_IN_FLUX, MF_COL_RANGE_FLUX);
            cpl_table_erase_column(rangespec, MF_COL_WEIGHT);
            cpl_table_name_column( rangespec, MF_COL_MOD_SCALE, MF_COL_RANGE_SCALE);
            cpl_table_erase_column(rangespec, MF_COL_MOD_FLUX);
            cpl_table_erase_column(rangespec, MF_COL_MOD_WEIGHT);
            cpl_table_erase_column(rangespec, MF_COL_DEV);

            /* Extract chip-related wavelength range from model spectrum */
            double limlam[2];

            limlam[0] = cpl_table_get(rangespec, MF_COL_IN_LAMBDA, 0,      NULL);
            limlam[1] = cpl_table_get(rangespec, MF_COL_IN_LAMBDA, nsel-1, NULL);

            limlam[0] = 1e4 / (1e4 / limlam[0] + MF_EXTRA_WN_COVERAGE);
            limlam[1] = 1e4 / (1e4 / limlam[1] - MF_EXTRA_WN_COVERAGE);

            cpl_table_or_selected_double( modspec, MF_COL_IN_LAMBDA, CPL_GREATER_THAN,     limlam[0]);
            cpl_table_and_selected_double(modspec, MF_COL_IN_LAMBDA, CPL_NOT_GREATER_THAN, limlam[1]);

            cpl_table *extmodspec = cpl_table_extract_selected(modspec);

            /* Delete model spectrum if no more required */
            if (!single_spectrum || (single_spectrum && j == nrange - 1)) {
                cpl_table_delete(modspec);
            }

            /* Get chip number for selected fit range */
            int chip = cpl_table_get(params->rangetab, MF_COL_CHIP, j, NULL);

            /* Modify wavelength grid of model spectrum by means of Chebyshev polynomials */
            mf_convolution_mod_wave_grid(extmodspec, params, fitpar, chip);

            /* Rebin model spectrum to wavelength grid of observed spectrum */
            mf_convolution_rebin(rangespec, MF_COL_IN_LAMBDA, MF_COL_IN_FLUX,    extmodspec, MF_COL_IN_LAMBDA, MF_COL_IN_FLUX );
            mf_convolution_rebin(rangespec, MF_COL_IN_LAMBDA, MF_COL_MOD_LAMBDA, extmodspec, MF_COL_IN_LAMBDA, MF_COL_LAMBDA_0);

            cpl_table_delete(extmodspec);

            /* Kernel spectrum convolution */
            if (kernel) mf_kernel_user(     rangespec, kernel, selrows);
            else        mf_kernel_synthetic(rangespec, params, fitpar );

            /* Add grey body in the case of a radiance spectrum and Adapt flux units */
            if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE && correct_spectrum) {
                mf_convolution_telback(rangespec, params, fitpar);
            }

            /* Modify continuum of model spectrum by means of polynomials */

            double flux0V[nsel];
            for (cpl_size i = 0; i < nsel; i++) {
                flux0V[i]=cpl_table_get(rangespec, MF_COL_IN_FLUX,     i, NULL);
            }
    
            cpl_table *rangespec_cpy =cpl_table_duplicate(rangespec);
            mf_convolution_mod_continuum(rangespec, params, fitpar, j + 1);

            /* Write resulting spectrum in output CPL table */
            cpl_table_duplicate_column(rangespec, "FLUX0", rangespec, MF_COL_IN_FLUX);

            for (cpl_size i = 0; i < nsel; i++) {

                cpl_size idx = cpl_array_get(selrows, i, NULL);

                cpl_table_set(spec, MF_COL_MOD_LAMBDA, idx, cpl_table_get(rangespec, MF_COL_MOD_LAMBDA,  i, NULL));
                cpl_table_set(spec, MF_COL_MOD_SCALE,  idx, cpl_table_get(rangespec, MF_COL_RANGE_SCALE, i, NULL));
                cpl_table_set(spec, MF_COL_MOD_FLUX,   idx, cpl_table_get(rangespec, MF_COL_IN_FLUX,     i, NULL));

                /* MNB*/
                /*cpl_table_dump_structure(rangespec,NULL);*/

                cpl_table_set(flux0, "FLUX0",           idx, cpl_table_get(rangespec, "FLUX0",     i, NULL));

            }

            cpl_msg_info(cpl_func,"RANGE %lld", j);

            int ncont  = 1 + params->config->fitting.fit_continuum.n;
            cpl_matrix * A   = cpl_matrix_new(nsel,ncont);
            cpl_matrix * RHS = cpl_matrix_new(nsel,1);
            cpl_matrix * SOL = cpl_matrix_new(ncont,1);
            cpl_matrix * W   = cpl_matrix_new(nsel,nsel);

            /* Shift zero point of wavelength scale to center of wavelength range */
            cpl_size nlam = cpl_table_get_nrow( rangespec);

            double limlam2[2] = { cpl_table_get( rangespec, MF_COL_IN_LAMBDA,  0,        NULL),
            cpl_table_get( rangespec, MF_COL_IN_LAMBDA,  nlam - 1, NULL) };

            double wmean = (limlam2[0] + limlam2[1]) / 2;

            cpl_msg_info(cpl_func,"MNB lam0=%f lam1=%f wmean=%f",limlam[0],limlam[1],wmean);
            cpl_size k,l;
            for( k=0;k<nsel;k++) {
                for (l=0;l<nsel;l++) {
                    cpl_matrix_set(W,k,l,0.0);
                }
            }
            for (k=0;k<nsel;k++)  {
                cpl_size idx2 = cpl_array_get(selrows, k, NULL);
                double val; /*,wt,wt1,wt2;*/
                   
                /*val = cpl_table_get(rangespec, MF_COL_IN_FLUX,k, NULL);*/
                val = cpl_table_get(spec, MF_COL_MOD_FLUX,idx2, NULL);
                cpl_matrix_set(RHS,k,0,val);
                /*
                val = cpl_table_get(spec, MF_COL_MOD_FLUX,idx2, NULL);
                */
                double wt1 = cpl_table_get(spec, MF_COL_WEIGHT,idx2, NULL);
                /* double wt2 = cpl_table_get(spec, MF_COL_MOD_WEIGHT,idx2, NULL);*/
                /*wt=sqrt(wt1*wt2);*/
                double wt=wt1;
                /*
                cpl_matrix_set(RHS,k,0,val);
                */
                cpl_matrix_set(W  ,k,k,wt);

            }


            for (l=0;l<ncont;l++) cpl_matrix_set(SOL,l,0,0.0);
                
            for (k=0;k<nsel;k++)  {
                double lam=cpl_table_get(rangespec, MF_COL_MOD_LAMBDA,  k, NULL);
                lam=lam-wmean;
                double val=flux0V[k];
                cpl_matrix_set(A,k,0,val);
                for (l=1;l<ncont;l++) {
                    val=val*lam;
                    cpl_matrix_set(A,k,l,val);
                }
            }

            cpl_matrix* A2   =cpl_matrix_product_create(W,A);
            cpl_matrix* RHS2 =cpl_matrix_product_create(W,RHS);

            SOL=cpl_matrix_solve_svd(A,RHS);
		
            gsl_matrix * gsl_M = gsl_matrix_alloc (nsel, ncont);
            gsl_vector * gsl_R = gsl_vector_alloc (nsel);
            gsl_vector * gsl_X = gsl_vector_alloc (ncont);

            for (int ig = 0; ig < nsel; ig++)
                gsl_vector_set (gsl_R, ig, cpl_matrix_get(RHS,ig,0));
 
            for (int ig = 0; ig < nsel; ig++)
                for (int jg = 0; jg < ncont; jg++)
                    gsl_matrix_set (gsl_M, ig, jg,  cpl_matrix_get(A,ig,jg));
 
            gsl_matrix * gsl_V = gsl_matrix_alloc (ncont,ncont);
            gsl_vector * gsl_S = gsl_vector_alloc (ncont);
            gsl_vector * gsl_Work = gsl_vector_alloc (ncont);
                
            gsl_linalg_SV_decomp(gsl_M,gsl_V,gsl_S,gsl_Work);
            gsl_linalg_SV_solve(gsl_M,gsl_V,gsl_S,gsl_R,gsl_X);
            /*gsl_linalg_SV_solve*/
            for (int jg = 0; jg < ncont; jg++) {
                cpl_msg_info(cpl_func,"GSL SVD, %d, %f %f",jg,gsl_vector_get (gsl_X,jg),gsl_vector_get(gsl_S,jg));
                cpl_msg_info(cpl_func,"CPL SVD, %d, %f   ",jg,cpl_matrix_get(SOL,jg,0));
            }
            gsl_matrix_free (gsl_M);
            gsl_vector_free (gsl_R);



            SOL=cpl_matrix_solve_svd(A2,RHS2);
            cpl_matrix_delete(A2);
            cpl_matrix_delete(RHS2);
                
            cpl_matrix_dump(SOL,NULL);


            int nmolec =     params->config->internal.molecules.n_molec;
            int nchip  =     params->config->internal.nchip;
            int nwlc   = 1 + params->config->fitting.fit_wavelenght.n;
            const double *par = cpl_array_get_data_double_const(fitpar);
            cpl_array *fitpar2 = cpl_array_duplicate(fitpar);

            int n0 = nmolec + (nwlc * nchip) + (ncont * (j + 0));

            /*mf_convolution_mod_continuum(rangespec, params, fitpar, j + 1);*/

            for  (int l1=0;l1<ncont;l1++) {
                double sval=cpl_matrix_get(SOL,l1,0);
                double fval=par[n0 + l1];
                cpl_msg_info(cpl_func,"RANGE %lld i=%d sol=%f, cof=%f cp %f",j,l1,sval,fval,cpl_array_get_double(fitpar2,n0+l1,NULL));
                    cpl_array_set_double (fitpar2,n0+l1,sval);
            }
		
            cpl_array_set_double (fitpar2,n0+4,0.0);
            /*
            cpl_array_set_double (fitpar2,n0+5,0.0);
            */
            cpl_matrix_delete(A);
            cpl_matrix_delete(W);
            cpl_matrix_delete(RHS);
            cpl_matrix_delete(SOL);


            if (last_call && 1==0) {
                cpl_msg_info(cpl_func,"LAST CALL");

                /* Modify continuum of model spectrum by means of polynomials */
                mf_convolution_mod_continuum(rangespec_cpy, params, fitpar2, j + 1);


                /* Write resulting spectrum in output CPL table */
                for (cpl_size i = 0; i < nsel; i++) {

                    cpl_size idx = cpl_array_get(selrows, i, NULL);

                    cpl_table_set(spec, MF_COL_MOD_LAMBDA, idx, cpl_table_get(rangespec_cpy, MF_COL_MOD_LAMBDA,  i, NULL));
                    cpl_table_set(spec, MF_COL_MOD_SCALE,  idx, cpl_table_get(rangespec_cpy, MF_COL_RANGE_SCALE, i, NULL));
                    cpl_table_set(spec, MF_COL_MOD_FLUX,   idx, cpl_table_get(rangespec_cpy, MF_COL_IN_FLUX,     i, NULL));

                    /* MNB*/
                    /*cpl_table_dump_structure(rangespec_cpy,NULL);*/
                }

                cpl_table_delete(rangespec);
                rangespec=rangespec_cpy;


            } else {
                cpl_table_delete(rangespec_cpy);
            };/* end if lats_call && 1==0)*/

            /* Cleanup */
            cpl_array_delete(selrows);
            cpl_table_delete(rangespec);
        }
    }; /*(cpl_size j = 0; j < nrange; j++*/
    
    for (cpl_size i = 0; i < m; i++) {

        int    range = cpl_table_get(spec, MF_COL_MOD_RANGE,  i, NULL);
        double mflux = cpl_table_get(spec, MF_COL_MOD_FLUX,   i, NULL);

        if (range == 0 || mflux < 0 || isnan(mflux) != 0) {
            cpl_table_set(spec, MF_COL_MOD_FLUX,   i, 0.);
            cpl_table_set(spec, MF_COL_MOD_WEIGHT, i, 0.);
        } else {
            cpl_table_set(spec, MF_COL_MOD_WEIGHT, i, 1.);
        }
    }

    /* Calculate weighted deviations between modelled and observed spectrum */
    cpl_table_add_columns(     spec, MF_COL_DEV, MF_COL_IN_FLUX  ); //cpl_table_add_columns(     spec, MF_COL_DEV, MF_COL_MOD_FLUX  );
    cpl_table_subtract_columns(spec, MF_COL_DEV, MF_COL_MOD_FLUX    ); //cpl_table_subtract_columns(spec, MF_COL_DEV, MF_COL_IN_FLUX   );
    cpl_table_multiply_columns(spec, MF_COL_DEV, MF_COL_WEIGHT    );
    cpl_table_multiply_columns(spec, MF_COL_DEV, MF_COL_MOD_WEIGHT);

    /* Print chi^2 (only for mf_model(...) not mf_calctrans...(...) */
    double chi2 = 0.;
    for (cpl_size i = 0; i < m; i++) {
        double dev = cpl_table_get(spec, MF_COL_DEV, i, NULL);
        chi2 += dev * dev;
    }


     /* MNB*/
     cpl_table_delete(flux0);

    if (last_call) cpl_msg_info(cpl_func, "(mf_convolutio) Last Iteration                     => Chi2: %12.6f",              chi2);
    else           cpl_msg_info(cpl_func, "(mf_convolutio) Loop Iteration (mpfit_calls = %4d) => Chi2: %12.6f", mpfit_calls, chi2);

    return CPL_ERROR_NONE;
}


/** @cond PRIVATE */

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_convolution_rebin(
    cpl_table                *outspec,
    const char               *outlam,
    const char               *outflux,
    const cpl_table          *inspec,
    const char               *inlam,
    const char               *influx)
{
    /*!
     * Rebins CPL table inspec with columns inlam and influx to wavelength
     * given by outlam in CPL table outspec. The resulting rebinned flux is
     * written to outflux in outspec. If outflux does not exist, it is
     * created. The routine conserves integral of flux.
     *
     * \b INPUT:
     * \param outspec  output spectrum with desired wavelength grid
     * \param outlam   output wavelength column name
     * \param outflux  output flux column name
     * \param inspec   input spectrum
     * \param inlam    input wavelength column name
     * \param influx   input flux column name
     *
     *
     * \b OUTPUT:
     * \param outspec  rebinned input spectrum
     *
     * \b ERRORS:
     * - none
     */

    /* Get number of data points */
    cpl_size n_in  = cpl_table_get_nrow(inspec );
    cpl_size n_out = cpl_table_get_nrow(outspec);

    /* Create flux column in output CPL table if it is not present */
    if (cpl_table_has_column(outspec, outflux) == 0) {
        cpl_table_new_column(outspec, outflux, CPL_TYPE_DOUBLE);
    }

    /* Fill flux column of output spectrum with 0. */
    cpl_table_fill_column_window(outspec, outflux, 0, n_out, 0.);

    /* No data points -> no rebinning */
    if (n_in <= 0 || n_out <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Get pointers to CPL table columns */
    const double *inlamv   = cpl_table_get_data_double_const(inspec,  inlam  );
    const double *influxv  = cpl_table_get_data_double_const(inspec,  influx );
    double       *outlamv  = cpl_table_get_data_double(      outspec, outlam );
    double       *outfluxv = cpl_table_get_data_double(      outspec, outflux);

    /* One input data point only */
    if (n_in == 1) {
        cpl_table_fill_column_window(outspec, outflux, 0, n_out, influxv[0]);
        return CPL_ERROR_NONE;
    }


    double   ilmin =  0.;
    double   ilmax =  0.;
    double   olmin =  0.;
    double   olmax =  0.;
    cpl_size jo    = -1;
    cpl_size j     =  0;
    for (cpl_size i = 0; i < n_out; i++) {

        /* Limits of wavelength bin in output spectrum */
        if (n_out == 1) {

          /* Full range of input spectrum for one output data point */
          olmin = 1.5 * inlamv[0]        - 0.5 * inlamv[1];
          olmax = 1.5 * inlamv[n_in - 1] - 0.5 * inlamv[n_in - 2];

        } else {

            if (i == 0)         olmin = 1.5 *  outlamv[i] - 0.5 * outlamv[i + 1];
            else                olmin = olmax;

            if (i == n_out - 1) olmax = 1.5 *  outlamv[i] - 0.5 * outlamv[i - 1];
            else                olmax = 0.5 * (outlamv[i] +       outlamv[i + 1]);
        }

        double dol = olmax - olmin;

        do {

            /* Limits of wavelength bin in input spectrum */
            if (j != jo) {

                if (j == 0)        ilmin = 1.5 *  inlamv[j] - 0.5 * inlamv[j + 1];
                else               ilmin = ilmax;

                if (j == n_in - 1) ilmax = 1.5 *  inlamv[j] - 0.5 * inlamv[j - 1];
                else               ilmax = 0.5 * (inlamv[j] +       inlamv[j + 1]);
            }

            /* Effective range of flux value -> weight */
            double dil;
            if (     ilmin <  olmin && ilmax <= olmax) dil = ilmax - olmin;
            else if (ilmin >= olmin && ilmax >  olmax) dil = olmax - ilmin;
            else if (ilmin <  olmin && ilmax >  olmax) dil = olmax - olmin;
            else                                       dil = ilmax - ilmin;

            /* Average flux of input spectrum in output bin */
            if (dil > 0) {
                double rdl = dil / dol;
                outfluxv[i] += influxv[j] * rdl;
            }

            j++;

        } while (ilmax <= olmax && j < n_in);

        j--;
        jo = j;
    }

    return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_convolution_mod_wave_grid(
    cpl_table                *spec,
    const mf_parameters      *params,
    const cpl_array          *fitpar,
    const int                chip)
{
     /*!
     * Modifies the wavelength scale of a model spectrum by means of
     * Chebyshev polynomials. The coefficients are provided by the fit
     * parameter vector (CPL array). The adaption of the wavelength grid is
     * carried out for the given chip only. The input wavelength grid is
     * saved and written into the output column MF_COL_LAMBDA_0.
     *
     * \b INPUT:
     * \param spec    CPL table with model spectrum
     * \param params  mf_parameters parameter structure
     * \param fitpar  CPL array with fit parameters
     * \param chip    chip number
     *
     * \b OUTPUT:
     * \param spec    model spectrum with modified wavelength grid
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    /* Find position of wavelength calibration parameters in fit parameter CPL array */
    int nwlc   = 1 + params->config->fitting.fit_wavelenght.n;
    int nmolec =     params->config->internal.molecules.n_molec;

    int n0 = nmolec + nwlc * (chip - 1);

    if (nwlc-1 < 1) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): wlc_n of mf_parameters *params < 1");
    }

    /* Get pointer to CPL array with fit parameters */
    const double *par = cpl_array_get_data_double_const(fitpar);

    /* Save input wavelength grid */
    cpl_table_duplicate_column(spec, MF_COL_LAMBDA_0, spec, MF_COL_IN_LAMBDA);

    /* Scale wavelengths in the way that the interval [-1,1] is covered */
    double limlam[2] = { cpl_table_get(params->chiptab, MF_COL_WL_MIN, chip-1, NULL),
                         cpl_table_get(params->chiptab, MF_COL_WL_MAX, chip-1, NULL) };

    double dlam  =  limlam[1] - limlam[0];
    double wmean = (limlam[0] + limlam[1]) / 2;

    cpl_table_duplicate_column(spec, MF_COL_SLAMBDA, spec, MF_COL_IN_LAMBDA);
    cpl_table_subtract_scalar( spec, MF_COL_SLAMBDA, wmean                 );
    cpl_table_divide_scalar(   spec, MF_COL_SLAMBDA, dlam / 2.             );

    /* Initialize wavelength column (all values = 0) */
    cpl_size nlam = cpl_table_get_nrow(spec);
    cpl_table_fill_column_window(spec, MF_COL_IN_LAMBDA, 0, nlam, 0.);

    /* Get pointers to CPL table columns */
    double *lam  = cpl_table_get_data_double(spec, MF_COL_IN_LAMBDA);
    double *slam = cpl_table_get_data_double(spec, MF_COL_SLAMBDA);

    /* Create array for Chebyshev polynomials */
    cpl_array *cheby = cpl_array_new(nwlc, CPL_TYPE_DOUBLE);
    double    *t     = cpl_array_get_data_double(cheby);

    /* Compute Chebyshev polynomials */
    t[0]  = 1;
    for (cpl_size i = 0; i < nlam; i++) {
        for (cpl_size j = 0; j < nwlc; j++) {

            if (     j == 1) t[j] =     slam[i];
            else if (j >  1) t[j] = 2 * slam[i] * t[j - 1] - t[j - 2];

            lam[i] += par[n0 + j] * t[j];
        }
    }

    /* Rescale wavelengths */
    cpl_table_multiply_scalar(spec, MF_COL_IN_LAMBDA, dlam / 2.);
    cpl_table_add_scalar(     spec, MF_COL_IN_LAMBDA, wmean    );

    /* Cleanup */
    cpl_array_delete(cheby);
    cpl_table_erase_column(spec, MF_COL_SLAMBDA);

    return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_convolution_telback(
    cpl_table                *spec,
    const mf_parameters      *params,
    const cpl_array          *fitpar)
{
    /*!
     * Computes thermal emission by telesope/instrument and adds it to the
     * input spectrum. The routine assumes a grey body depending on emissivity
     * (fit parameter telback) and ambient temperature (from mf_parameters
     * structure).
     *
     * \note The grey body emission is divided by 1 - emissivity in order to
     * have the same flux level as the sky emission. This approach is correct
     * for flux-calibrated spectra, where the reflection of light by the
     * telescope mirror has been corrected, which is wrong for the telescope
     * emission. For this reason, the apparent telescope emission becomes
     * higher than the true one.
     *
     * \b INPUT:
     * \param spec    input spectrum (CPL table)
     * \param params  mf_parameters parameter structure
     * \param fitpar  CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param spec    input spectrum + telescope/instrument emission
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    /* Constants and unit conversions */
    double c1 = 2 * CPL_PHYS_H * pow(CPL_PHYS_C, 2);
    double c2 = CPL_PHYS_H * CPL_PHYS_C / CPL_PHYS_K;
    double c3 = MF_CONVOLUTION_LAM_UNIT / (CPL_PHYS_H * CPL_PHYS_C * MF_CONVOLUTION_SR_IN_ARCSEC2);
    double c4 = c1 * c3 / pow(MF_CONVOLUTION_LAM_UNIT, 4);
    double c5 = c2 / MF_CONVOLUTION_LAM_UNIT;

    /* Get primary mirror temperature */
    double mirror_temperature = params->config->ambient.mirror_temperature.value + MF_CONV_T_2_K;
    if (mirror_temperature < 0.) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): m1temp of mf_parameters *params < 0 K");
    }

    /* Get emissivity estimate from CPL array of fit parameters */
    cpl_size ntel = cpl_array_get_size(fitpar) - 1;
    double eps = cpl_array_get(fitpar, ntel, NULL);
    if (eps < 0 || eps >= 1.) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): emissivity of cpl_array *fitpar < 0 or >= 1");
    }

    /* Get pointers to CPL table columns */
    double *lam  = cpl_table_get_data_double(spec, MF_COL_IN_LAMBDA);
    double *flux = cpl_table_get_data_double(spec, MF_COL_IN_FLUX  );

    /* Compute and add grey body emission to spectrum */
    cpl_size nlam = cpl_table_get_nrow(spec);
    for (cpl_size i = 0; i < nlam; i++) {
        flux[i] += (eps / (1 - eps)) * c4 / (pow(lam[i], 4) * (expm1(c5 / (lam[i] * mirror_temperature))));
    }

    /* Conversion of fluxes from phot/(s*m2*mum*as2) (emission spectrum only)
     * to flux unit of observed spectrum.
     * Meaning of the flag flux_unit of the mf_parameters parameter structure:
     * 0: phot/(s*m^2*mum*as^2) [no conversion]
     * 1: W/(m^2*mum*as^2)
     * 2: erg/(s*cm^2*A*as^2)
     * 3: mJy/as^2
     * For other units, the conversion factor has to be considered as constant
     * term of the continuum fit parameters.
     */

    /* Get flux unit flag */
    int flux_unit = params->config->fitting.flux_unit;

    /* Convert fluxes */
    if (       flux_unit == MF_PARAMETERS_FLUX_UNIT_1) {

        /* phot/sm2mum -> W/m2mum */
        cpl_table_divide_columns(  spec, MF_COL_IN_FLUX, MF_COL_IN_LAMBDA);
        cpl_table_multiply_scalar( spec, MF_COL_IN_FLUX,        CPL_PHYS_C * CPL_PHYS_H / MF_CONVOLUTION_LAM_UNIT);

    } else if (flux_unit == MF_PARAMETERS_FLUX_UNIT_2) {

        /* phot/sm2mum -> erg/scm2A */
        cpl_table_divide_columns(  spec, MF_COL_IN_FLUX, MF_COL_IN_LAMBDA);
        cpl_table_multiply_scalar( spec, MF_COL_IN_FLUX,  0.1 * CPL_PHYS_C * CPL_PHYS_H / MF_CONVOLUTION_LAM_UNIT);

    } else if (flux_unit == MF_PARAMETERS_FLUX_UNIT_3) {

        /* phot/sm2mum -> mJy */
        cpl_table_multiply_columns(spec, MF_COL_IN_FLUX, MF_COL_IN_LAMBDA);
        cpl_table_multiply_scalar( spec, MF_COL_IN_FLUX, 1e35 *              CPL_PHYS_H * MF_CONVOLUTION_LAM_UNIT );
    }

    return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_convolution_mod_continuum(
    cpl_table                *spec,
    const mf_parameters      *params,
    const cpl_array          *fitpar,
    const int                range)
{
    /*!
     * Modifies continuum of a model spectrum by means of polynomials.
     * The coefficients are provided by the fit parameter vector (CPL array).
     * The continuum correction is carried out for the given fit range only.
     *
     * \b INPUT:
     * \param spec    CPL table with model spectrum
     * \param params  mf_parameters parameter structure
     * \param fitpar  CPL array with fit parameters
     * \param range   range number
     *
     * \b OUTPUT:
     * \param spec    model spectrum with modified continuum
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    /* Find position of wavelength calibration parameters in fit parameter CPL array */
    int nmolec =     params->config->internal.molecules.n_molec;
    int nchip  =     params->config->internal.nchip;
    int ncont  = 1 + params->config->fitting.fit_continuum.n;
    int nwlc   = 1 + params->config->fitting.fit_wavelenght.n;

    int n0 = nmolec + (nwlc * nchip) + (ncont * (range - 1));
    if (ncont - 1 < 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): cont_n of mf_parameters *params < 0");
    }

    /* Get pointer to CPL array with fit parameters */
    const double *par = cpl_array_get_data_double_const(fitpar);

    /* Shift zero point of wavelength scale to center of wavelength range */
    cpl_size nlam = cpl_table_get_nrow( spec);

    double limlam[2] = { cpl_table_get( spec, MF_COL_IN_LAMBDA,  0,        NULL),
                         cpl_table_get( spec, MF_COL_IN_LAMBDA,  nlam - 1, NULL) };

    double wmean = (limlam[0] + limlam[1]) / 2;

    cpl_table_duplicate_column(spec, MF_COL_SLAMBDA, spec,   MF_COL_IN_LAMBDA);
    cpl_table_subtract_scalar( spec, MF_COL_SLAMBDA, wmean);

    /* Get pointers to CPL table columns */
    double *slam = cpl_table_get_data_double(spec, MF_COL_SLAMBDA    );
    double *scal = cpl_table_get_data_double(spec, MF_COL_RANGE_SCALE);

/* MNB FROM HERE */

                cpl_matrix * A   = cpl_matrix_new(nlam,ncont);
                cpl_matrix * RHS = cpl_matrix_new(nlam,1);
                cpl_matrix * SOL = cpl_matrix_new(ncont,1);
/* MNB TO HERE */

    /* Compute continuum correction polynomials */
    for (cpl_size i = 0; i < nlam; i++) {
                 /* MNB FROm HERE */
                 double val = 0.0; /*cpl_table_get(spec, MF_COL_MOD_FLUX,i, NULL);*/
                 cpl_matrix_set(RHS,i,0,val);
		 /* MNB TO HERE */
 
        double fac = 0.;
        for (cpl_size j = 0; j < ncont; j++) {
            fac += par[n0 + j] * pow(slam[i], j);
	    
                 /* MNB FROm HERE */
		 double flx0 =cpl_table_get(spec, MF_COL_IN_FLUX,i,NULL);
                 cpl_matrix_set(A,i,j,fac*flx0);
		 /* MNB TO HERE */
        }

        if (fac <= 0) scal[i] = 1.;
        else          scal[i] = fac;
    }

    /* Remove temporary table column */
    cpl_table_erase_column(spec, MF_COL_SLAMBDA);

    /* Scale continuum (but save unscaled for later analysis MNB) */
    cpl_table_multiply_columns(spec, MF_COL_IN_FLUX, MF_COL_RANGE_SCALE);
    
/* MNB FROM HERE */
        for (cpl_size j = 0; j < ncont; j++) {
	    cpl_msg_info(cpl_func,"STORED Coeff %lld = %f",j,par[n0+j]);

        }               
            cpl_msg_info(cpl_func,"BEG ---- MNB HERE 2 -------->");
	        SOL=cpl_matrix_solve_svd(A,RHS);
            cpl_matrix_dump(SOL,NULL);
            cpl_msg_info(cpl_func,"END ---- MNB HERE 2 -------->");

		/*cpl_matrix_dump(SOL,NULL);*/
                cpl_matrix_delete(A);
                cpl_matrix_delete(RHS);
                cpl_matrix_delete(SOL);
               
/* MNB TO HERE */

    return CPL_ERROR_NONE;
}

/** @endcond */


/**@}*/
