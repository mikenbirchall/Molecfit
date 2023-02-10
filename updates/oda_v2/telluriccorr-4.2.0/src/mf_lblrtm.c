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

#include "mf_io.h"

#include "mf_lblrtm.h"

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

/* Execute the lblrtm binary in all of the wavenumbers in a concrete range */
static cpl_error_code mf_lblrtm_range_execution(
    mf_io_lblrtm_config      *config_lblrtm,
    mf_parameters            *params,
    const int                ismolcol,
    const cpl_size           run_execution,
    int                      *lblrtm_calls,
    cpl_table                *atm_profile,
    cpl_table                **spec_out,
    cpl_vector               *mol_abuns,
    cpl_error_code           *range_status);

/*  */
static cpl_error_code mf_lblrtm_molecular_correct_atm(
    cpl_table                *atm_profile,
    cpl_table                **atm_profile_out,
    const mf_parameters      *params,
    const cpl_array          *fitpar);

/* Make output spectrum : Combined and convolved */
static cpl_table * mf_lblrtm_range_combined_and_convolved(
    const char               *lblrtm_out_filename,
    const double             pixel_res,
    const double             wn1,
    const double             wn2,
    const char               *w_dir_range, const int range, cpl_vector* mol_abuns, cpl_boolean USE_HYBRID);

/*  */
static cpl_error_code mf_lblrtm_create_wavelength_grid(
    cpl_table                *outspec,
    const double             limlam[2],
    const double             resol);

/* Spectrum linear interpolation */
static cpl_error_code mf_lblrtm_linear_interpolate_spectrum(
    double                   *lambda,
    double                   *flux,
    const cpl_size           lim_min,
    const cpl_size           lim_max);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_lblrtm       .
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create the user parameter structure (with default values).
 *
 * @param .                  .
 * @param                    .
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
mf_io_lblrtm_config * mf_lblrtm_config_create(void)
{
    mf_io_lblrtm_config *config = cpl_malloc(sizeof(mf_io_lblrtm_config));

    config->icntnm    = MF_LBLRTM_ICNTNM_INIT;
    config->iaersl    = MF_LBLRTM_IAERSL_INIT;
    config->mpts      = MF_LBLRTM_MPTS_INIT;
    config->npts      = MF_LBLRTM_NPTS_INIT;
    config->v[0]      = MF_LBLRTM_V1_INIT;
    config->v[1]      = MF_LBLRTM_V2_INIT;
    config->sample    = MF_LBLRTM_SAMPLE_INIT;
    config->alfal0    = MF_LBLRTM_ALFAL0_INIT;
    config->avmass    = MF_LBLRTM_AVMASS_INIT;
    config->dptmin    = MF_LBLRTM_DPTMIN_INIT;
    config->dptfac    = MF_LBLRTM_DPTFAC_INIT;
    config->tbound    = MF_LBLRTM_TBOUND_INIT;
    config->sremis[0] = MF_LBLRTM_SREMIS1_INIT;
    config->sremis[1] = MF_LBLRTM_SREMIS2_INIT;
    config->sremis[2] = MF_LBLRTM_SREMIS3_INIT;
    config->srrefl[0] = MF_LBLRTM_SRREFL1_INIT;
    config->srrefl[1] = MF_LBLRTM_SRREFL2_INIT;
    config->srrefl[2] = MF_LBLRTM_SRREFL3_INIT;
    config->model     = MF_LBLRTM_MODEL_INIT;
    config->itype     = MF_LBLRTM_ITYPE_INIT;
    config->nozero    = MF_LBLRTM_NOZERO_INIT;
    config->noprnt    = MF_LBLRTM_NOPRNT_INIT;
    config->ipunch    = MF_LBLRTM_IPUNCH_INIT;
    config->re        = MF_LBLRTM_RE_INIT;
    config->hspace    = MF_LBLRTM_HSPACE_INIT;
    config->ref_lat   = MF_PARAMETERS_LATITUDE_VALUE_INIT; //MF_LBLRTM_REF_LAT_INIT;
    config->h[0]      = MF_PARAMETERS_ELEVATION_VALUE_INIT/1000.0; // MF_LBLRTM_H1_INIT;
    config->h[1]      = MF_LBLRTM_H2_INIT;
    config->range     = MF_LBLRTM_RANGE_INIT;
    config->beta      = MF_LBLRTM_BETA_INIT;
    config->len       = MF_LBLRTM_LEN_INIT;
    config->hobs      = MF_LBLRTM_HOBS_INIT;
    config->avtrat    = MF_LBLRTM_AVTRAT_INIT;
    config->tdiff[0]  = MF_LBLRTM_TDIFF1_INIT;
    config->tdiff[1]  = MF_LBLRTM_TDIFF2_INIT;
    config->altd[0]   = MF_LBLRTM_ALTD1_INIT;
    config->altd[1]   = MF_LBLRTM_ALTD2_INIT;
    config->delv      = MF_LBLRTM_DELV_INIT;

    return config;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Deallocate the user parameter structure.
 *
 * @param config             .
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
void mf_lblrtm_config_delete(
    mf_io_lblrtm_config      *config)
{
  if (config) {

      cpl_free(config);
  }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Check the user parameter structure. For check the user modifications.
 *
 * @param config             .
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
cpl_error_code mf_lblrtm_config_check(
    mf_io_lblrtm_config      *config)
{
  cpl_ensure(config, CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);

  // TODO: Add check for all variables and values
  cpl_error_code err = CPL_ERROR_NONE;

  if (!err && (   config->icntnm < MF_LBLRTM_ICNTNM_MIN
               || config->icntnm > MF_LBLRTM_ICNTNM_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_ICNTNM);
  }

  if (!err && (   config->iaersl < MF_LBLRTM_IAERSL_MIN
               || config->iaersl > MF_LBLRTM_IAERSL_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_IAERSL);
  }

  if (!err && (   config->sample < MF_LBLRTM_SAMPLE_MIN
               || config->sample > MF_LBLRTM_SAMPLE_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_SAMPLE);
  }

  if (!err && (   config->model < MF_LBLRTM_MODEL_MIN
               || config->model > MF_LBLRTM_MODEL_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_MODEL);
  }

  if (!err && (   config->itype < MF_LBLRTM_ITYPE_MIN
               || config->itype > MF_LBLRTM_ITYPE_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_ITYPE);
  }

  if (!err && (   config->nozero < MF_LBLRTM_NOZERO_MIN
               || config->nozero > MF_LBLRTM_NOZERO_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_NOZERO);
  }

  if (!err && (   config->noprnt < MF_LBLRTM_NOPRNT_MIN
               || config->noprnt > MF_LBLRTM_NOPRNT_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_NOPRNT);
  }

  if (!err && (   config->ipunch < MF_LBLRTM_IPUNCH_MIN
               || config->ipunch > MF_LBLRTM_IPUNCH_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_IPUNCH);
  }

  if (!err && (   config->ref_lat < MF_PARAMETERS_LATITUDE_VALUE_MIN // MF_LBLRTM_REF_LAT_MIN
               || config->ref_lat > MF_PARAMETERS_LATITUDE_VALUE_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_LATITUDE_VALUE);
  }

  if (!err && (   config->len < MF_LBLRTM_LEN_MIN
               || config->len > MF_LBLRTM_LEN_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_LEN);
  }


  return err;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief Execute LBLRTM
 *
 * @param run_execution      .
 * @param params             .
 * @param last_call          .
 * @param spec_out           .
 * @param range_status       .
 * @param reffitpar          .
 * @param mpfit_calls        .
 * @param lblrtm_calls       .
 * @param config_lblrtm      .
 * @param atm_profile        .
 * @param atm_profile_out    .
 * @param fitpar             .
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
cpl_error_code mf_lblrtm(
    const cpl_size           run_execution,
    mf_parameters            *params,
    cpl_boolean              last_call,
    cpl_table                **spec_out,
    cpl_error_code           *range_status,
    cpl_array                *reffitpar,
    int                      mpfit_calls,
    int                      *lblrtm_calls,
    mf_io_lblrtm_config      *config_lblrtm,
    const cpl_table          *atm_profile,
    cpl_table                **atm_profile_out,
    cpl_array                *fitpar)
{
    /*** First Part of the call to LBLRTM ***/

    /* Check whether a relative molecular column density was changed (only for mf_model(...) not mf_calctrans...(...) */
    int ismolcol = 1;
    if (mpfit_calls != 1 || !last_call) {

        /* Check number of array values */
         cpl_size n    = cpl_array_get_size(fitpar);
         cpl_size nref = cpl_array_get_size(reffitpar);
         if (n == 0 && n == nref) {
             ismolcol = 1;
             return CPL_ERROR_NONE;
         } else if (n != nref) {
             ismolcol = -1;
             return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                          "Invalid object structure: cpl_array *fitpar != cpl_array *reffitpar (size)");
         }

         /* Get number of molecules in driver parameter structure */
         int nmolec = params->config->internal.molecules.n_molec;

         /* Check modification of molecular column density parameter and copy fitpar to reffitpar */
         if (mpfit_calls == 1) cpl_msg_debug(cpl_func, "(mf_lblrtm    ) FIT ITERATION: %4d ... Getting first approximation",                mpfit_calls        );
         else                  cpl_msg_debug(cpl_func, "(mf_lblrtm    ) FIT ITERATION: %4d ... Checking tolerance [MF_TOLERANCE =  %5.2e]", mpfit_calls, MF_TOL);

         cpl_boolean all_in_tolerance = CPL_TRUE;

         ismolcol = 0;
         for (cpl_size i = 0; i < n; i++) {

             double par    = cpl_array_get(fitpar,    i, NULL);
             double refpar = cpl_array_get(reffitpar, i, NULL);

             double reldev = fabs(par - refpar);

             if (     refpar == 0. && par == 0.) reldev  = 0.;
             else if (refpar == 0.             ) reldev /= par;
             else                                reldev /= refpar;

             if (mpfit_calls != 1 && reldev > MF_TOL) {
                 cpl_msg_debug(cpl_func, "(mf_lblrtm    ) [ELEM = %3lld, VAL = %14.5f, REF = %14.5f] => |VAL-REF|/VAL = %5.2e > %5.2e (NOT FITTING)",
                               i, par, refpar, reldev, MF_TOL);
                 if (all_in_tolerance) all_in_tolerance = CPL_FALSE;
             }

             if (i < nmolec && reldev > MF_TOL) ismolcol = 1;

             /* Check if reffitpar have enough elements */
             if (cpl_array_get_size(reffitpar) > i) {
                 cpl_array_set(reffitpar, i, par);
             }
         }

         if (mpfit_calls != 1 && all_in_tolerance) cpl_msg_info(cpl_func, "(mf_lblrtm    ) ALL FITTING !");
    }

    /* Run LBLRTM in any case if first or last function call */
    if (mpfit_calls == 1 || last_call) ismolcol = 1;

    /* Modify profiles of atmospheric molecules if necessary */
    cpl_table *atm_profile_local = NULL;
    if (ismolcol == 1) {
        atm_profile_local = cpl_table_duplicate(atm_profile);
        mf_lblrtm_molecular_correct_atm(atm_profile_local, atm_profile_out, params, fitpar);
    }

    /* Execute the LBLRTM calls */
    cpl_vector *mol_abuns = mf_io_molecule_abundancies(params,fitpar);
    cpl_error_code err = mf_lblrtm_range_execution(config_lblrtm, params,
                                                   ismolcol,
                                                   run_execution, lblrtm_calls,
                                                   atm_profile_local,
                                                   spec_out, mol_abuns, range_status);

    cpl_vector_delete(mol_abuns);
    /* Cleanup */
    if (atm_profile_local) cpl_table_delete(atm_profile_local);

    return err;
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
static cpl_error_code mf_lblrtm_molecular_correct_atm(
    cpl_table                *atm_profile,
    cpl_table                **atm_profile_out,
    const mf_parameters      *params,
    const cpl_array          *fitpar)
{
    /*!
     * \callgraph
     *
     * Multiplies profiles of molecules listed in the mf_parameters structure by
     * correction factors given by the CPL array of fit parameters.
     * The resulting CPL table of atmospheric profiles is written into an
     * ASCII output file (parameter MF_PARAMETERS_OUTPUT_NAME in mf_parameters structure +
     * MF_FILE_FIT_ATM_ASCII) which is required as input file for the LBLRTM code.
     *
     * \b INPUT:
     * \param modprof  CPL table with atmospheric profiles
     * \param params   mf_parameters parameter structure
     * \param fitpar   CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param modprof  atmospheric profiles with corrected molecular
     *                 frequencies
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    /* Get number of columns with molecular profiles */
    cpl_size nmolcol = cpl_table_get_ncol(atm_profile) - 3;

    /* Get pointer to column names */
    cpl_array *colarray = cpl_table_get_column_names(atm_profile);
    char      **colname = cpl_array_get_data_string(colarray);

    /* Get number of molecules in driver parameter structure */
    int nmolec = params->config->internal.molecules.n_molec;

    /* Get pointer to names of molecules in driver parameter structure */
    char **mol = cpl_table_get_data_string(params->molectab, MF_COL_LIST_MOLECULES);

    /* Get pointer to CPL array with fit parameters */
    const double *par = cpl_array_get_data_double_const(fitpar);

    /* Modify column densities of molecules */
    cpl_error_code err = CPL_ERROR_NONE;
    for (cpl_size i = 0; i < nmolec; i++) {
        cpl_boolean ext = CPL_FALSE;
        for (cpl_size j = 0; j < nmolcol; j++) {
           if (strcmp(colname[j + 3], mol[i]) == 0) {
                ext = CPL_TRUE;
                cpl_table_multiply_scalar(atm_profile, colname[j + 3], par[i]);
            }
        }

        /* Molecule of the driver file is not present in list of CPL table column names */
        if (!ext) err = cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                              "Invalid object structure: cpl_table *modprof (column %s not found)",
                                              mol[i]);
    }

    /* Create output table. If not the first call, rewrite it. Otherwise, we end up with a memory leak. */
    if (*atm_profile_out) cpl_table_delete(*atm_profile_out);
    *atm_profile_out = cpl_table_duplicate(atm_profile);

    /* Hack for the first value of the pressure */
    cpl_table_set_double(*atm_profile_out, MF_COL_ATM_PRE, 0, cpl_table_get_double(atm_profile, MF_COL_ATM_PRE, 0, NULL) + 0.1);

    /* Cleanup */
    cpl_array_delete(colarray);

    return err;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Execute the LBLRTM binary in all of the wavenumbers in a concrete range.
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
static cpl_error_code mf_lblrtm_range_execution(
    mf_io_lblrtm_config      *config_lblrtm,
    mf_parameters            *params,
    const int                ismolcol,
    const cpl_size           run_execution,
    int                      *lblrtm_calls,
    cpl_table                *atm_profile,
    cpl_table                **spec_out,
    cpl_vector               *mol_abuns,
    cpl_error_code           *range_status)
{

    /* Local parameters */
    cpl_error_code err               = CPL_ERROR_NONE;
    cpl_boolean    no_single_spec    = !params->config->internal.single_spectrum;
    int            nrange            =  params->config->internal.n_range;

    char           **run_w_dir_range = cpl_calloc(nrange, sizeof(char *));
    double         *wn_end           = cpl_calloc(nrange, sizeof(double));
    double         *pixel_res        = cpl_calloc(nrange, sizeof(double));
    cpl_boolean    *execute_range    = cpl_calloc(nrange, sizeof(cpl_boolean));

    cpl_boolean USE_HYBRID=CPL_TRUE;
    //cpl_boolean USE_HYBRID=CPL_FALSE;

    if (params->config->internal.single_spectrum) {
        cpl_msg_info(cpl_func, "(mf_lblrtm    ) Compute single spectrum with mf_lblrtm(...)");
    }

    /* Create folders */
    for (cpl_size range = 0; range < nrange; range++) {

        /* Compute LBLRTM spectrum. Get range (number of LBLRTM spectra) */
        run_w_dir_range[range] = cpl_sprintf("%s/run_%lld_%s_%lld_lblrtm_call_%lld",
                                             params->config->internal.tmp_folder,
                                             run_execution,
                                             MF_AER_WDIR_LBLRTM_RANGE_PATH,
                                             range + 1,
                                             *lblrtm_calls + range + 1);

        /* Get wavenumber range in [cm-1] */
        wn_end[   range] = cpl_table_get(params->rangetab, MF_COL_WN_END,  range, NULL)             ;
        pixel_res[range] = cpl_table_get(params->rangetab, MF_COL_PIX_RES, range, NULL) * MF_OVERSAMPLING_FACTOR;

        if (range == 0 || no_single_spec) {

            /* Skip empty ranges */
            if (wn_end[range] == 0.) {

                spec_out[range]     = NULL;
                range_status[range] = CPL_ERROR_NONE;

            } else if (ismolcol == 1) {

                if (pixel_res[range] == 0.) {

                    spec_out[range]     = NULL;
                    range_status[range] = CPL_ERROR_NONE;

                } else {

                    /* Create working folders */
                    if (mf_io_mkdir(run_w_dir_range[range]) != CPL_ERROR_NONE) {
                        err += cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                                        "Could not create working directory");
                    }
                }
            }
        }
    }


    if (!err) {

        /* Get name LBLRTM output : Emission/Transmission */
        char *lblrtm_out_filename;
        cpl_boolean spec_emission = params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE;
        if (spec_emission) lblrtm_out_filename = cpl_sprintf("%s", MF_AER_TAPE27_FILE);
        else               lblrtm_out_filename = cpl_sprintf("%s", MF_AER_TAPE28_FILE);

        double   *wn1                 = cpl_calloc(nrange, sizeof(double    ));
        double   *wn2                 = cpl_calloc(nrange, sizeof(double    ));
        cpl_size *num_wavenumber      = cpl_calloc(nrange, sizeof(cpl_size  ));
        cpl_size **name_wavenumber    = cpl_calloc(nrange, sizeof(cpl_size *));
        char     ***folder_wavenumber = cpl_calloc(nrange, sizeof(char *    ));

        /* Execute ranges and update the number of lblrtm function calls */
        for (cpl_size range = 0; range < nrange; range++) {

            /* Get wavenumber range */
            wn1[range] = cpl_table_get(params->rangetab, MF_COL_WN_START, range, NULL);
            wn2[range] = cpl_table_get(params->rangetab, MF_COL_WN_END,   range, NULL);

            /* Prepare the execution files */
            if ( (range == 0 || no_single_spec)  && wn_end[range] != 0. && ismolcol == 1 && pixel_res[range] != 0. ) {

                cpl_size wn_jmax = 0;
                double   min_wn  = wn2[range];
                double   max_wn  = min_wn;
                for (cpl_size j = 0; min_wn > wn1[range]; j++) {

                    double delta = 0.14 * pow(max_wn, 1.1);
                    if (delta > MF_LBLRTM_DELTA_MAX) delta = MF_LBLRTM_DELTA_MAX;
                    if (delta < MF_LBLRTM_DELTA_MIN) delta = MF_LBLRTM_DELTA_MIN;

                    min_wn = CPL_MAX(max_wn - delta, wn1[range]);
                    max_wn = min_wn;

                    wn_jmax = j;
                }

                /* Calculation of frequency for refractive geometry. It need to convert wavelenght to wavenumber and calculate */
                double minc_config = MF_CONV_K_LAM / config_lblrtm->v[0];
                double maxc_config = MF_CONV_K_LAM / config_lblrtm->v[1];

                /* Check wavelenght range */
                if (minc_config < 0 || maxc_config < 0 || fabs(maxc_config - minc_config) > MF_LBLRTM_DELTA_ABS) {
                    err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                                "LBLRTM variable delta too big (|%g - %g| = %g) => |maxc - minc| > %g cm-1",
                                                maxc_config, minc_config, fabs(maxc_config - minc_config), MF_LBLRTM_DELTA_ABS);
                }

                double     vbar        = (minc_config + maxc_config) / 2.;
                double     angle       = 90. - params->config->ambient.telescope_angle.value;
                const char *lbl_molecs = params->config->internal.molecules.lbl_molecs;

                /* Ranges */
                double         minc;
                double         maxc;

                num_wavenumber[   range] = 0;
                name_wavenumber[  range] = NULL;
                folder_wavenumber[range] = NULL;

                for (cpl_size i = 0; i < 2; i++) {

                    /* 1er run : Count the number of valid wavenumbers */
                    /* 2nd run : Execute the configuration             */

                    if (i != 0) {
                        name_wavenumber[  range] = cpl_calloc(num_wavenumber[range], sizeof(cpl_size));
                        folder_wavenumber[range] = cpl_calloc(num_wavenumber[range], sizeof(char *  ));
                    }

                    cpl_size j            = -1;
                    double   min_wn_local = wn2[range];
                    double   max_wn_local = min_wn_local;

                    for (cpl_size wavenumber = wn_jmax; min_wn_local > wn1[range] && !err; wavenumber--) {

                        double delta = MF_LBLRTM_DELTA_FACTOR * pow(max_wn_local, MF_LBLRTM_DELTA_TOL);
                        if (delta > MF_LBLRTM_DELTA_MAX) delta = MF_LBLRTM_DELTA_MAX;
                        if (delta < MF_LBLRTM_DELTA_MIN) delta = MF_LBLRTM_DELTA_MIN;

                        min_wn_local = CPL_MAX(max_wn_local - delta, wn1[range]);

                        if (wavenumber == wn_jmax) maxc = max_wn_local;
                        else                       maxc = max_wn_local + MF_EXTRA_WN_COVERAGE;

                        if (wavenumber == 0      ) minc = min_wn_local;
                        else                       minc = min_wn_local - MF_EXTRA_WN_COVERAGE;

                        /* Update for the next loop */
                        max_wn_local = min_wn_local;

                        if (i == 0) {

                            /* Increase number of waves */
                            num_wavenumber[range]++;

                        } else {

                            /* Save the concrete number of wave */
                            j++;
                            name_wavenumber[  range][j] = wavenumber;
                            folder_wavenumber[range][j] = cpl_sprintf("%s/wv_number_%lld", run_w_dir_range[range], wavenumber);

                            /* Check inputs parameters */
                            if (minc < 0 || minc >= maxc || fabs(maxc - minc) > MF_LBLRTM_DELTA_ABS) {

                                err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                                             "(|%g - %g| = %g) => |maxc - minc| > %g cm-1",
                                                             maxc, minc, fabs(maxc - minc), MF_LBLRTM_DELTA_ABS);
                            } else {

                                err = mf_io_mkdir(folder_wavenumber[range][j]);
                                if (!err) {

                                      /* Create a symbolic link to TAPE3 (output of LNFL) and write the TAPE5 input file by LNFL execution */
                                      const char *tape3 = cpl_table_get_string(params->rangetab, MF_COL_LNFL, range);
                                      err = mf_io_write_lblrtm_configuration(folder_wavenumber[range][j], tape3,
                                                                             minc, maxc,
                                                                             vbar, angle, spec_emission, lbl_molecs,
                                                                             config_lblrtm,
                                                                             atm_profile);
                                }
                            }
                        }
                    }
                }

                execute_range[range] = (!err) ? CPL_TRUE : CPL_FALSE;

            } else {

                execute_range[range] = CPL_FALSE;
            }
        }

        if (!err) {

            cpl_msg_info(cpl_func, "(mf_lblrtm    ) Executing external calls to the LBLRTM binary ...");

            char *lblrtm_syscall;
            if (params->config->inputs.silent_external_bins) {
                lblrtm_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LBLRTM_SILENT);
            } else {
                lblrtm_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LBLRTM);
            }

            int    num_calls      = 0;
            double runtime_lblrtm = 0.;
            #ifdef _USE_OPENMP
                /* OMP parallel version */
                /*if (params->config->inputs.omp_num_threads > 1) {
                    cpl_msg_info(cpl_func, "(mf_lblrtm   ) OPENMP in parallel by wavelength ranges");
                }*/
                #pragma omp parallel for shared(lblrtm_syscall, execute_range, num_wavenumber, folder_wavenumber) reduction(+ : err, num_calls, runtime_lblrtm)
            #endif
            for (cpl_size range = 0; range < nrange; range++) {
                if (execute_range[range]) {
                    for (cpl_size wavenumber = 0; wavenumber < num_wavenumber[range]; wavenumber++) {
                        double runtime;
                        cpl_matrix* odTable =mf_io_oda_get_tableDB(range);
                        if (!USE_HYBRID || (USE_HYBRID && odTable==NULL)) {
                        #ifdef _USE_OPENMP
                            err += mf_io_systemOpenMP(lblrtm_syscall, folder_wavenumber[range][wavenumber], &runtime);
                        #else
                            err += mf_io_system      (lblrtm_syscall, folder_wavenumber[range][wavenumber], &runtime);
                        #endif
                        runtime_lblrtm += runtime;
                        }
                    }
                    num_calls++;
                }
            }
            *lblrtm_calls                  += num_calls;
            params->timers.time_bin_lblrtm += runtime_lblrtm;

            mf_io_sync();

            cpl_free(lblrtm_syscall);


            if (!err) {

                for (cpl_size range = 0; range < nrange && !err; range++) {

                      if (execute_range[range]) {

                          /* Rename and move output file */
                          for (cpl_size wavenumber = 0; wavenumber < num_wavenumber[range]; wavenumber++) {

                              char *old_output_name = cpl_sprintf("%s/%s",      folder_wavenumber[range][wavenumber], lblrtm_out_filename                                        );
                              char *new_output_name = cpl_sprintf("%s/%s_%lld", run_w_dir_range[range],               lblrtm_out_filename, name_wavenumber[range][wavenumber] + 1);
                                cpl_matrix* odTable =mf_io_oda_get_tableDB(range);
                                if (!USE_HYBRID || (USE_HYBRID && odTable==NULL)) {
                                    err = mf_io_mv(old_output_name, new_output_name);
                                }
									
			      if (err!=0) exit(1);

                              cpl_free(old_output_name);
                              cpl_free(new_output_name);
                          }


                          if (err != CPL_ERROR_NONE) {

                              cpl_error_set_message(cpl_func, err,
                                                    "LBLRTM failed for wavelength interval: %5.2f-%5.2f [µm]",
                                                    MF_CONV_K_LAM / wn1[range], MF_CONV_K_LAM / wn2[range]);
                          } else {

                              /* Get pixel resolution */
                              double local_pixel_res = CPL_MIN(1e6, pixel_res[range]);

                              /* Create spectrums combine the outputs for all the wavenumbers */
                              spec_out[range] = mf_lblrtm_range_combined_and_convolved(lblrtm_out_filename, local_pixel_res, wn1[range], wn2[range], run_w_dir_range[range],range,mol_abuns, USE_HYBRID);

                              /* Save output error */
                              err = (spec_out[range]) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_OUTPUT;
                              range_status[range] = err;
                          }
                      }
                }
            }
        }

        /* Cleanup */
        for (cpl_size range = 0; range < nrange; range++) {

            if (name_wavenumber[range]) cpl_free(name_wavenumber[range]);

            if (folder_wavenumber[range]) {
                for (cpl_size wavenumber = 0; wavenumber < num_wavenumber[range]; wavenumber++){
                    cpl_free(folder_wavenumber[range][wavenumber]);
                }
                cpl_free(folder_wavenumber[range]);
            }
        }
        cpl_free(num_wavenumber);
        cpl_free(name_wavenumber);
        cpl_free(folder_wavenumber);
        cpl_free(wn1);
        cpl_free(wn2);
        cpl_free(lblrtm_out_filename);
    }

    /* Cleanup */
    for (cpl_size range = 0; range < nrange; range++) {
        cpl_free(run_w_dir_range[range]);
    }
    cpl_free(run_w_dir_range);
    cpl_free(wn_end         );
    cpl_free(pixel_res      );
    cpl_free(execute_range  );

    return err;
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
static cpl_table * mf_lblrtm_range_combined_and_convolved(
    const char               *lblrtm_out_filename,
    const double             pixel_res,
    const double             wn1,
    const double             wn2,
    const char               *w_dir_range,const int range, cpl_vector* mol_abuns, cpl_boolean USE_HYBRID)
{
    /* Merges and rebins a set of LBLRTM radiance ("*R.*") or transmission ("*T.*") spectra */
    double limlam[2] = { MF_CONV_K_LAM / wn2,
                         MF_CONV_K_LAM / wn1 };

    /* Initialize wavelength grid */
    cpl_table *spec = cpl_table_new(0);
    cpl_table_new_column(spec, MF_COL_IN_LAMBDA, CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, MF_COL_IN_FLUX,   CPL_TYPE_DOUBLE);

    /* Create grid */
    cpl_error_code err = mf_lblrtm_create_wavelength_grid(spec, limlam, pixel_res);
    if (err != CPL_ERROR_NONE) {
        cpl_error_set_message(cpl_func, err,
                              "Error in mf_lib_createlibspec");
        cpl_table_delete(spec);
        return NULL;
    }

    /*** Rebin set of LBLRTM spectra (wrapper output) in wavelength units [mu m] (variable step size possible) **/

    /* Find number of files and wavenumber ranges. */
    cpl_array *klim_all;
    if (USE_HYBRID) {
        klim_all=mf_io_klim_from_odatable(range);
        if (klim_all==NULL) klim_all = mf_io_find_klim(w_dir_range, lblrtm_out_filename);
    } else {
        klim_all = mf_io_find_klim(w_dir_range, lblrtm_out_filename);
    }
    if (!klim_all) {
        err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                    "Invalid input parameter(s): klim_all NULL!");
    } else {

        /* Rebin spectra */
        cpl_size       size = cpl_array_get_size(klim_all);
        for (cpl_size wavenumber = 0; wavenumber < size - 1 && !err; wavenumber++) {

            double klim[2] = { cpl_array_get(klim_all, wavenumber,     NULL),
                               cpl_array_get(klim_all, wavenumber + 1, NULL)  };

            /* Check wavenumber range */
            if (klim[1] <= klim[0]) {
                err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                             "Invalid input parameter(s): klim[1] <= klim[0] (wavenumber limits)");
            } else {

                /* Convert wavenumbers to wavelengths */
                double llim[2] = { MF_CONV_K_LAM / klim[1],
                                   MF_CONV_K_LAM / klim[0] };

                /* Check number of data points in spectrum */
                int nrow = cpl_table_get_nrow(spec);
                if (nrow == 0) {
                    err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                                "No data: cpl_table *spec");
                }

                /* Get pointers to CPL table columns */
                double *lamv  = cpl_table_get_data_double(spec, MF_COL_IN_LAMBDA);
                double *fluxv = cpl_table_get_data_double(spec, MF_COL_IN_FLUX  );

                /* Compare wavelength ranges */
                if (llim[1] < lamv[0] || llim[0] > lamv[nrow - 1]) {

                    /* Desired wavelength range not covered by the input model spectrum */
                    err = CPL_ERROR_NONE;

                } else {

                    /* Output variable of Rebin -> mf_io_read_lblrtm_and_update_spec(...) */
                    cpl_boolean usampl =  CPL_FALSE;
                    int         jmin   = -1;
                    int         jmax   =  0;

                    /* Rebin spectrum */
                    char *spectrum_filename = cpl_sprintf("%s/%s_%lld", w_dir_range, lblrtm_out_filename, wavenumber + 1);
                    cpl_bivector* bvec;
                    bvec=mf_io_mergeODTables(range,mol_abuns,spectrum_filename);
                    if (bvec!=NULL) {
                        if (USE_HYBRID) {
                            cpl_msg_info(cpl_func,"MNBXX:BVEC was not null and not deleting to NULL");
                        } else {
                            cpl_bivector_delete(bvec);
                            bvec=NULL;
                            cpl_msg_info(cpl_func,"MNBXX:BVEC was not null so deleted and reset to NULL");                        }
                    }
                    if (bvec==NULL) bvec = mf_io_read_lblrtm_spec(spectrum_filename);
                    err = mf_io_read_lblrtm_and_update_spec(nrow, lamv, fluxv, bvec, llim, &usampl, &jmin, &jmax);
                    cpl_free(spectrum_filename);
                    cpl_bivector_delete(bvec);

                    /* Interpolate "empty" bins by means of the closest bins that contain data points (valid flux values) */
                    if (usampl) mf_lblrtm_linear_interpolate_spectrum(lamv, fluxv, jmin, jmax);
                }
            }
        }
    }

    /* Cleanup */
    cpl_array_delete(klim_all);

    /* Check error and return */
    if (!err) {
        return spec;
    } else {
      cpl_table_delete(spec);
      return NULL;
    }
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
static cpl_error_code mf_lblrtm_create_wavelength_grid(
    cpl_table                *outspec,
    const double             limlam[2],
    const double             resol)
{
    /*!
     * Initialisation of a wavelength grid that is characterised by constant
     * resolution lambda / delta lambda
     *
     * \b INPUT:
     * \param limlam  lower and upper limit of wavelength range
     * \param resol   resolution lambda / delta lambda
     *
     * \b OUTPUT:
     * \param outspec  CPL table that provides wavelengths and fluxes (= 0.)
     *                 of a spectrum
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    /* Test input parameters */
    if (limlam[1] <= limlam[0] || resol <= 0.) {
        /* Return spectrum with zero data points */
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid input parameter(s): limlam[1] <= limlam[0] || resol <= 0. (wavelength grid)");
    }

    /* Get number of data points and resize CPL table */
    cpl_size nrow = ceil(log(limlam[1] / limlam[0]) / log1p(1 / resol));
    cpl_table_set_size(outspec, nrow);

    /*** Create wavelength grid ***/

    /* Set first position */
    double lam = limlam[0];
    cpl_table_set(outspec, MF_COL_IN_LAMBDA, 0, lam);

    /* Set next positions */
    for (cpl_size i = 1; i < nrow; i++) {
        lam *= (1 + 1 / resol);
        cpl_table_set(outspec, MF_COL_IN_LAMBDA, i, lam);
    }

    /* Set flux = 0 */
    cpl_table_fill_column_window(outspec, MF_COL_IN_FLUX, 0, nrow, 0.);

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
static cpl_error_code mf_lblrtm_linear_interpolate_spectrum(
    double                   *lambda,
    double                   *flux,
    const cpl_size           lim_min,
    const cpl_size           lim_max)
{
  /*!
   * Linear interpolation for "empty" bins in spectra (marked by -HUGE_VAL).
   *
   * \b INPUT:
   * \param lim_min  pixel min limit of spectrum
   * \param lim_max  pixel max limit of spectrum
   *
   * \b OUTPUT:
   * \param spec    spectrum corrected by interpolation
   *
   * \b ERRORS:
   * - No data
   */

  /* Perform interpolation/extrapolation procedure */
  int            jlim[2] = {-1, -1};
  cpl_boolean    ok      = CPL_TRUE;
  cpl_error_code status  = CPL_ERROR_NONE;

  for (cpl_size j = lim_min; j <= lim_max; j++) {

      if (flux[j] != -HUGE_VAL) {

          if (ok == CPL_TRUE) {
              jlim[0] = j;
          } else {
              jlim[1] = j;
              if (jlim[0] == -1) {
                  /* Constant value for extrapolation of lower margin */
                  for (cpl_size i = lim_min; i <= jlim[1] - 1; i++) {
                      flux[i] = flux[jlim[1]];
                  }
              } else {
                  /* Linear interpolation */
                  double m = (flux[jlim[1]] - flux[jlim[0]]) / (lambda[jlim[1]] - lambda[jlim[0]]);
                  for (cpl_size i = jlim[0] + 1; i <= jlim[1] - 1; i++) {
                      flux[i] = m * (lambda[i] - lambda[jlim[0]]) + flux[jlim[0]];
                  }
              }
              jlim[0] = jlim[1];
              ok = CPL_TRUE;
          }

      } else {

          ok = CPL_FALSE;

          if (j == lim_max) {
              if (jlim[0] == -1) {
                  status = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                                 "No data: No valid flux value");
              } else {
                  /* Constant value for extrapolation of upper margin */
                  for (cpl_size i = jlim[0] + 1; i <= lim_max; i++) {
                      flux[i] = flux[jlim[0]];
                  }
              }
          }

      }

  }

  return status;
}

/** @endcond */


/**@}*/
