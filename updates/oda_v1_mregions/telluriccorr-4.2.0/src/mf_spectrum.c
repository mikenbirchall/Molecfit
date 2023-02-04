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

#include "mf_spectrum.h"

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
static cpl_error_code mf_spectrum_include_ranges(
    mf_parameters            *params,
    cpl_table                *spec,
    const cpl_table          *ranges);

/*  */
static cpl_error_code mf_spectrum_exclude_ranges(
    cpl_table                *spec,
    const cpl_table          *ranges,
    const char               wp);

/*  */
static cpl_error_code mf_spectrum_air_to_vacuum(
    const char               *frame_type,
    const double             obs_erf_rv,
    cpl_table                *spectrum);

/*  */
static cpl_error_code mf_spectrum_fill_parameter_table(
    mf_parameters            *params,
    cpl_table                *spectrum);

/*  */
static cpl_error_code mf_spectrum_replace_coefficients(
    mf_parameters            *params,
    cpl_table                *rangetab,
    cpl_table                *chiptab);

/*  */
static cpl_error_code mf_spectrum_set_weights(
    const mf_parameters      *params,
    cpl_table                *spec);

/*  */
static cpl_error_code mf_spectrum_calculate_wave_range(
    mf_parameters            *params,
    cpl_table                *spec);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_spectrum       .
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create input spectrum from call telluriccorr
 *
 * @param spectrum           .
 * @param column_lambda      .
 * @param column_flux        .
 * @param column_dflux       .
 * @param column_mask        .
 * @param telluric_corr      .
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
cpl_table * mf_spectrum_create_input_to_correct(
    const cpl_table          *spectrum,
    const char               *column_lambda,
    const char               *column_flux,
    const char               *column_dflux,
    const char               *column_mask,
    const cpl_vector         *telluric_corr)
{
  /* Check inputs */
  cpl_error_ensure(spectrum && telluric_corr, CPL_ERROR_NULL_INPUT, return NULL, "NULL input");

  /* Check number  of data points (check zeroth extension only) */
  cpl_size n_row_spectrum = cpl_table_get_nrow(spectrum);
  cpl_error_ensure(n_row_spectrum > 0, CPL_ERROR_ILLEGAL_INPUT, return NULL, "Not rows in spectrum");

  /* Create spectrum table */
  cpl_table *telluriccorr_data = cpl_table_new(n_row_spectrum);
  cpl_table_new_column(telluriccorr_data, MF_COL_CHIP,              CPL_TYPE_INT   );
  cpl_table_new_column(telluriccorr_data, MF_COL_IN_LAMBDA,         CPL_TYPE_DOUBLE);
  cpl_table_new_column(telluriccorr_data, MF_COL_IN_FLUX,           CPL_TYPE_DOUBLE);
  cpl_table_new_column(telluriccorr_data, MF_COL_IN_DFLUX,          CPL_TYPE_DOUBLE);
  cpl_table_new_column(telluriccorr_data, MF_COL_IN_MASK,           CPL_TYPE_INT   );
  cpl_table_new_column(telluriccorr_data, MF_COL_OUT_FLUX,          CPL_TYPE_DOUBLE);
  cpl_table_new_column(telluriccorr_data, MF_COL_OUT_DFLUX,         CPL_TYPE_DOUBLE);
  cpl_table_new_column(telluriccorr_data, MF_COL_OUT_MASK,          CPL_TYPE_INT   );
  cpl_table_new_column(telluriccorr_data, MF_COL_OUT_TELLURIC_CORR, CPL_TYPE_DOUBLE);

  const char  *column_lambda_str = strcmp(column_lambda, MF_PARAMETERS_COLUMN_DFLUX_NULL) ? column_lambda : MF_PARAMETERS_COLUMN_LAMBDA_DEFAULT;
  const char  *column_flux_str   = strcmp(column_flux,   MF_PARAMETERS_COLUMN_FLUX_NULL ) ? column_flux   : MF_PARAMETERS_COLUMN_FLUX_DEFAULT;
  cpl_boolean exist_column_dflux = strcmp(column_dflux,  MF_PARAMETERS_COLUMN_DFLUX_NULL);
  cpl_boolean exist_column_mask  = strcmp(column_mask,   MF_PARAMETERS_COLUMN_MASK_NULL );

  cpl_error_code err = cpl_error_get_code();

  /* Copy data */
  for (cpl_size i = 0; i < n_row_spectrum && !err; i++) {

      int    chip     =                      cpl_table_get(spectrum, MF_COL_CHIP,       i, NULL);
      double lambda   =                      cpl_table_get(spectrum, column_lambda_str, i, NULL);
      double flux     =                      cpl_table_get(spectrum, column_flux_str,   i, NULL);
      double dflux    = exist_column_dflux ? cpl_table_get(spectrum, column_dflux,      i, NULL) : 0.;             /* Default : Non-error      */
      double mask     = exist_column_mask  ? cpl_table_get(spectrum, column_mask,       i, NULL) : CPL_BINARY_0;   /* Default : Valid          */

      double corr     = telluric_corr      ? cpl_vector_get(telluric_corr, i)                    : 0.;             /* Default : Non-correction */

      cpl_table_set(telluriccorr_data, MF_COL_CHIP,              i, chip  );
      cpl_table_set(telluriccorr_data, MF_COL_IN_LAMBDA,         i, lambda);
      cpl_table_set(telluriccorr_data, MF_COL_IN_FLUX,           i, flux  );
      cpl_table_set(telluriccorr_data, MF_COL_IN_DFLUX,          i, dflux );
      cpl_table_set(telluriccorr_data, MF_COL_IN_MASK,           i, mask  );
      cpl_table_set(telluriccorr_data, MF_COL_OUT_FLUX,          i, flux  );
      cpl_table_set(telluriccorr_data, MF_COL_OUT_DFLUX,         i, dflux );
      cpl_table_set(telluriccorr_data, MF_COL_OUT_MASK,          i, mask  );
      cpl_table_set(telluriccorr_data, MF_COL_OUT_TELLURIC_CORR, i, corr  );

      err = cpl_error_get_code();
  }

  if (!err) {
      return telluriccorr_data;
  } else {
      cpl_table_delete(telluriccorr_data);

      return NULL;
  }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create input spectrum from call telluriccorr
 *
 * @param config             Configuration parameters for create the internal cpl_table
 * @param spectrum           Input spectrum data table
 *
 * @return cpl_table         data if everything is OK or NULL in other case.
 *
 * @note Puts the read data in a CPL table consisting of the columns MF_COL_LAMBDA, MF_COL_FLUX, and MF_COL_WEIGHT.
 *       The names of the required FITS table columns are provided by the mf_parameters_config structure.
 *       The presence of flux error is optional.
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_table * mf_spectrum_create(
    mf_parameters_config     *config,
    const cpl_table          *spectrum)
{
    /*!
     * Reads a FITS file with tabulated spectroscopic data (wavelength, flux,
     * flux error, mask) that was created by mf_spectrum_preptable and puts
     * the read data in a CPL table consisting of the columns MF_COL_LAMBDA,
     * MF_COL_FLUX, and MF_COL_WEIGHT. The names of the required FITS table columns are
     * provided by the mf_parameters parameter structure. The presence of flux error
     * is optional. For skipping such a column the name has to be 'NULL'. The
     * original input file could also miss a mask column (also indicated by
     * 'NULL'). However, mf_spectrum_preptable makes sure that a suitable
     * mask column indicated by the given column name or mf_DEFMASKCOL (in
     * the case of 'NULL') + '_I' is present. The input file can have an
     * arbitrary number of extensions with different parts of the spectrum (to
     * be defined in the telluriccorr driver file). A weight of zero is taken if a
     * wavelength is excluded by the mask.
     *
     */
     
    cpl_error_code err = CPL_ERROR_NONE;

    /*** COLUMNS ***/
    int ncolmin = 5;

    char *column_lambda = NULL;
    if (strcmp(config->inputs.column_lambda, MF_PARAMETERS_COLUMN_LAMBDA_NULL) == 0) {
        column_lambda = cpl_sprintf("%s", MF_PARAMETERS_COLUMN_LAMBDA_DEFAULT);
    } else {
        column_lambda = cpl_sprintf("%s", config->inputs.column_lambda);
    }

    char *column_flux = NULL;
    if (strcmp(config->inputs.column_flux, MF_PARAMETERS_COLUMN_FLUX_NULL) == 0) {
        column_flux = cpl_sprintf("%s", MF_PARAMETERS_COLUMN_FLUX_DEFAULT);
    } else {
        column_flux = cpl_sprintf("%s", config->inputs.column_flux);
    }

    char *column_dflux = cpl_sprintf("%s", config->inputs.column_dflux);
    cpl_boolean exist_column_dflux = CPL_TRUE;
    if (strcmp(column_dflux, MF_PARAMETERS_COLUMN_DFLUX_NULL) == 0) {
        ncolmin--;
        exist_column_dflux = CPL_FALSE;
    }

    char *column_mask  = cpl_sprintf("%s", config->inputs.column_mask);
    char *column_imask = NULL;
    cpl_boolean exist_column_mask = CPL_TRUE;
    if (strcmp(column_mask, MF_PARAMETERS_COLUMN_MASK_NULL) == 0) {
        ncolmin--;
        exist_column_mask = CPL_FALSE;
        column_imask = cpl_sprintf("%s%s", MF_PARAMETERS_COLUMN_MASK_DEFAULT, "_I");
    } else {
        column_imask = cpl_sprintf("%s%s", column_mask,                       "_I");
    }

    /* Check data for each chip */
    cpl_boolean isnanum    = CPL_FALSE;
    cpl_boolean isoutrange = CPL_FALSE;
    cpl_boolean isnomask   = CPL_FALSE;
    cpl_boolean ismask0    = CPL_FALSE;
    cpl_boolean is0        = CPL_FALSE;
    int         nmask0     = 0;

    /* Check existence of expected columns */
    err = CPL_ERROR_NONE;
    if (   cpl_table_has_column(spectrum, column_lambda)  != 1
        || cpl_table_has_column(spectrum, column_flux  )  != 1) {
        err = CPL_ERROR_INCOMPATIBLE_INPUT;
    } else if (exist_column_dflux){
        if (cpl_table_has_column(spectrum, column_dflux ) != 1) err = CPL_ERROR_INCOMPATIBLE_INPUT;
    } else if (exist_column_mask) {
        if (cpl_table_has_column(spectrum, column_mask  ) != 1) err = CPL_ERROR_INCOMPATIBLE_INPUT;
    }

    if (err != CPL_ERROR_NONE) {

        cpl_free(column_lambda);
        cpl_free(column_flux  );
        cpl_free(column_dflux );
        cpl_free(column_mask  );
        cpl_free(column_imask );

        cpl_msg_info(cpl_func, "Column search for wavelength = %s", config->inputs.column_lambda);
        cpl_msg_info(cpl_func, "Column search for flux       = %s", config->inputs.column_flux);
        cpl_msg_info(cpl_func, "Column search for dflux      = %s", config->inputs.column_dflux);
        cpl_msg_info(cpl_func, "Column search for mask       = %s", config->inputs.column_mask);

        cpl_error_set_message(cpl_func, err,
                                     "Invalid spectrum column structure");

        return NULL;
    }

    /* Copy input table */
    cpl_table *spec_adapt = cpl_table_duplicate(spectrum);

    if (cpl_table_has_column(spec_adapt, column_imask) == 1) {

        char *column_renamed = cpl_sprintf("%s_orig", column_imask);
        if (cpl_table_has_column(spec_adapt, column_renamed) != 1) {
           cpl_msg_info(cpl_func, "(mf_io        ) Name of internal integer mask already used: Rename %s in %s", column_imask, column_renamed);
           cpl_table_name_column(spec_adapt, column_imask, column_renamed);
        } else {
           cpl_msg_info(cpl_func, "(mf_io        ) Use of reserved mask column names: Erase %s, keep %s", column_imask, column_renamed);
           cpl_table_erase_column(spec_adapt, column_imask);
        }
        cpl_free(column_renamed);
    }
    cpl_table_new_column(spec_adapt, column_imask, CPL_TYPE_INT);

    /* Check for nan values and negative values */
    cpl_boolean isnanflux  = CPL_FALSE;
    cpl_boolean isnandflux = CPL_FALSE;

    /* Correct possible error values */
    cpl_size nrow = cpl_table_get_nrow(spec_adapt);
    for (cpl_size i = 0; i < nrow; i++) {

        double flux = cpl_table_get(spec_adapt, column_flux, i, NULL);

        if (isnan(flux) != 0) {
            cpl_table_set(spec_adapt, column_flux, i, 0.);
            isnanflux = CPL_TRUE;
            isnanum   = CPL_TRUE;
        } else {
            isnanflux = CPL_FALSE;
        }

        if (exist_column_dflux) {

            double dflux = cpl_table_get(spec_adapt, column_dflux, i, NULL);

            if (dflux <= 0 || isnan(dflux) != 0) {
                cpl_table_set(spec_adapt, column_dflux, i, 0.);
                isnandflux = CPL_TRUE;
                isoutrange = CPL_TRUE;
            } else {
                isnandflux = CPL_FALSE;
            }
        }

        double mask = exist_column_mask ? cpl_table_get(spec_adapt, column_mask, i, NULL) : 1.;

        if (isnanflux || isnandflux) {

            cpl_table_set(spec_adapt, column_imask, i, -1);
            if (mask != 0 && mask != 1) {
                isnomask = CPL_TRUE;
            } else if (mask == 0) {
                is0 = CPL_TRUE;
                nmask0++;
            }

        } else {

            if (mask == 0) {
                cpl_table_set(spec_adapt, column_imask, i, 0);
                is0 = CPL_TRUE;
                nmask0++;
            } else {
                cpl_table_set(spec_adapt, column_imask, i, 1);
                if (mask != 1) {
                    isnomask = CPL_TRUE;
                }
            }
        }
    }

    /* Return if mask cannot be interpreted */
    if (isnomask && !is0) {

        cpl_free(column_lambda);
        cpl_free(column_flux  );
        cpl_free(column_dflux );
        cpl_free(column_mask  );
        cpl_free(column_imask );

        cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                     "Invalid object structure: mf_data (mf_data_table) *tabdat (all mask value(s) != 0 or 1)");

        cpl_table_delete(spec_adapt);
        return NULL;
    }

    /* Correct integer mask values if required (e.g. reverse definition) */
    for (cpl_size i = 0; i < nrow; i++) {

        int imask = cpl_table_get(spec_adapt, column_imask, i, NULL);

        if (isnomask) {

            if (imask == 0) {
                cpl_table_set(spec_adapt, column_imask, i, 1);
            } else {
                cpl_table_set(spec_adapt, column_imask, i, 0);
            }

        } else {

            if (imask == -1) {
                cpl_table_set(spec_adapt, column_imask, i, 0);
            } else if (nmask0 == nrow) {
                cpl_table_set(spec_adapt, column_imask, i, 1);
                ismask0 = CPL_TRUE;
            }
        }
    }

    /* Print info message in the case of bad fluxes, errors, or mask values */
    if (isnanum   ) cpl_msg_info(cpl_func, "(mf_io        ) Input data: flux(es) = 'nan'        -> set mask = 0");
    if (isoutrange) cpl_msg_info(cpl_func, "(mf_io        ) Input data: error(s) <= 0 or 'nan'  -> set mask = 0");
    if (isnomask  ) cpl_msg_info(cpl_func, "(mf_io        ) Input data: mask value(s) != 0 or 1 -> reverse definition (0 -> 1; != 0 -> 0)");
    if (ismask0   ) cpl_msg_info(cpl_func, "(mf_io        ) Input data: all mask values = 0     -> reverse definition (0 -> 1)");


    /* Get wavelength unit conversion factor */
    double wlgtomicron = config->inputs.wlg_to_micron;


    /* Create output table columns */
    cpl_table *spec_out = cpl_table_new(0);
    cpl_table_new_column(spec_out, MF_COL_CHIP,       CPL_TYPE_INT   );
    cpl_table_new_column(spec_out, MF_COL_IN_LAMBDA,  CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec_out, MF_COL_IN_FLUX,    CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec_out, MF_COL_WEIGHT,     CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec_out, MF_COL_MOD_RANGE,  CPL_TYPE_INT   );
    cpl_table_new_column(spec_out, MF_COL_MOD_LAMBDA, CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec_out, MF_COL_MOD_SCALE,  CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec_out, MF_COL_MOD_FLUX,   CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec_out, MF_COL_MOD_WEIGHT, CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec_out, MF_COL_DEV,        CPL_TYPE_DOUBLE);


    /* Read selected extensions of FITS file, compute weights and fill output CPL table */
    double min_lam = 1e200;
    double max_lam = 0.;

    /* Get column labels */
    cpl_array *colnames = cpl_table_get_column_names(spec_adapt);
    cpl_size  ncol      = cpl_array_get_size(colnames);

    /* Check existence of columns */
    cpl_size check   = 0;
    int      coln[4] = {0, 0, 0, 0};
    for (cpl_size i = 0; i < ncol; i++) {

        const char *colname = cpl_array_get_string(colnames, i);

        if (     strcmp(colname, column_lambda) == 0) { coln[0] = i; check++; }
        else if (strcmp(colname, column_flux  ) == 0) { coln[1] = i; check++; }
        else if (strcmp(colname, column_dflux ) == 0) { coln[2] = i; check++; }
        else if (strcmp(colname, column_imask ) == 0) { coln[3] = i; check++; }
        else if (strcmp(colname, column_mask  ) == 0) { check++;              }
    }

    if (check < ncolmin) {
        cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                              "Missing cpl table column(s)");
        cpl_table_delete(spec_out);
        cpl_array_delete(colnames);
        cpl_table_delete(spec_adapt);
        cpl_free(column_lambda);
        cpl_free(column_flux  );
        cpl_free(column_dflux );
        cpl_free(column_mask  );
        cpl_free(column_imask );
        return NULL;
    }

    /* Resize output table */
    cpl_table_set_size(spec_out, nrow);

    /* Transfer wavelength and flux/transmission and compute weights */
    cpl_size    chip_max    = 1;
    cpl_boolean column_chip = cpl_table_has_column(spec_adapt, MF_COL_CHIP);
    for (cpl_size i = 0; i < nrow; i++) {

        cpl_size chip = column_chip ? cpl_table_get(spec_adapt, MF_COL_CHIP, i, NULL) : 1;

        if (chip > chip_max) chip_max = chip;

        double lam  = cpl_table_get(spec_adapt, cpl_array_get_string(colnames, coln[0]), i, NULL) * wlgtomicron;
        double flux = cpl_table_get(spec_adapt, cpl_array_get_string(colnames, coln[1]), i, NULL);

        cpl_table_set(spec_out, MF_COL_CHIP,      i, chip);
        cpl_table_set(spec_out, MF_COL_IN_LAMBDA, i, lam );
        cpl_table_set(spec_out, MF_COL_IN_FLUX,   i, flux);

        min_lam = CPL_MIN(lam, min_lam);
        max_lam = CPL_MAX(lam, max_lam);

        double dflux;
        if (exist_column_dflux) dflux = cpl_table_get(spec_adapt, cpl_array_get_string(colnames, coln[2]), i, NULL);
        else                    dflux = 1.;

        int mask = cpl_table_get(spec_adapt, cpl_array_get_string(colnames, coln[3]), i, NULL);

        if (dflux <= 0. || mask == 0) cpl_table_set(spec_out, MF_COL_WEIGHT, i, 0.        );
        else                          cpl_table_set(spec_out, MF_COL_WEIGHT, i, 1. / dflux);
    }

    /* Delete temporary CPL objects */
    cpl_array_delete(colnames);
    cpl_table_delete(spec_adapt);


    /* Update nchips */
    config->internal.nchip = chip_max;

    /**********************************/

    /* Cleanup */
    cpl_free(column_lambda);
    cpl_free(column_flux  );
    cpl_free(column_dflux );
    cpl_free(column_mask  );
    cpl_free(column_imask );

    /* Check wavelength range, warn only (PIPE-5418) */
    if (min_lam < MF_WAVELENGTH_MIN_MICRONS || max_lam > MF_WAVELENGTH_MAX_MICRONS) {
        cpl_msg_warning(cpl_func,
                        "(mf_spectrum  ) Wavelength range is %g to %g micron (with wlgtomicron=%g). Is this really intended?",
                        min_lam, max_lam, wlgtomicron);
    } else if (min_lam <= 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                                     "Negative wavelength range, wlgtomicron=%g)",
                                     wlgtomicron);
        cpl_table_delete(spec_out);
        return NULL;
    }

    return spec_out;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Check ranges cpl_table
 *
 * @param spectrum           .
 * @param table              .
 * @param tag                .
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
cpl_error_code mf_spectrum_ranges_check(
    const cpl_table          *spectrum,
    const cpl_table          *table,
    const char               *tag)
{


  cpl_ensure(table, CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);

  cpl_error_ensure(   cpl_table_has_column(table, MF_COL_WAVE_RANGE_LOWER) != 0
                   && cpl_table_has_column(table, MF_COL_WAVE_RANGE_UPPER) != 0
                   && cpl_table_get_nrow(table)                             > 0,
                   CPL_ERROR_INCOMPATIBLE_INPUT, return CPL_ERROR_INCOMPATIBLE_INPUT,
                   "cpl_table *%s does not have the correct columns [%s, %s]",
                   tag, MF_COL_WAVE_RANGE_LOWER, MF_COL_WAVE_RANGE_UPPER);

  if (       !strcmp(tag, MF_INPUT_WAVE_INCLUDE)) {


 
      double lambda_min = cpl_table_get_column_min(spectrum, MF_COL_IN_LAMBDA);
      double lambda_max = cpl_table_get_column_max(spectrum, MF_COL_IN_LAMBDA);
      
      if (lambda_min>lambda_max) return CPL_ERROR_INCOMPATIBLE_INPUT;
      
     /* Anti-Sanity Check Requirement requested by Alain. i.e. remove this check */
/*

      double range_min  = cpl_table_get_column_min(table, MF_COL_WAVE_RANGE_LOWER);
      double range_max  = cpl_table_get_column_max(table, MF_COL_WAVE_RANGE_UPPER);

      cpl_error_ensure(   lambda_min <= range_min
                       && lambda_max >= range_max,
                       CPL_ERROR_INCOMPATIBLE_INPUT, return CPL_ERROR_INCOMPATIBLE_INPUT,
                       "The real spectrum cpl_table wavelengths range are %s[min=%g, max=%g] but the cpl_table *%s ranges are outside this range [%s(min)=%g, %s(max)=%g]",
                       MF_COL_IN_LAMBDA, lambda_min, lambda_max, tag, MF_COL_WAVE_RANGE_LOWER, range_min, MF_COL_WAVE_RANGE_UPPER, range_max);
*/
  } else if (!strcmp(tag, MF_INPUT_WAVE_EXCLUDE)) {

      /* No check need : If outside of wave range nothing happen */

  } else if (!strcmp(tag, MF_INPUT_PIXEL_EXCLUDE)) {

      /* No check need : If outside of pixel range nothing happen */
  }

  return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create the telluriccorr internal cpl_table spectrum structure.
 *
 * @param params             mf_parameters structure
 * @param spectrum           Data   spectrum
 * @param inc_wranges        Include wavelengths
 * @param exc_wranges        Exclude wavelengths
 * @param exc_pranges        Exclude pixels
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @note Reads a FITS table that was created by PREPTABLE and puts the read data
 *         in a CPL table consisting of the columns MF_COL_LAMBDA, MF_COL_FLUX, and
 *         MF_COL_WEIGHT. An input FITS file can have an arbitrary number of extensions
 *         with different parts of the spectrum (to be defined in the telluriccorr
 *         driver file). If the hierarchical ESO keywords are given in the header
 *         of the provided FITS file, the required parameters describing the
 *         telescope site and observing conditions are taken from the header.
 *         Their values substitute those given in the telluriccorr driver file.
 *         Moreover, files consisting of wavelength ranges to be included or
 *         wavelength or pixel ranges to be excluded are optionally read. These
 *         data determine the model range as given by the MF_SPEC_COL_MOD_RANGE column and the
 *         range table of the mf_parameters parameter structure.
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_table * mf_spectrum_ranges_apply(
    mf_parameters            *params,
    const cpl_table          *spectrum,
    const cpl_table          *inc_wranges,
    const cpl_table          *exc_wranges,
    const cpl_table          *exc_pranges)
{
    /* Check inputs */
    cpl_ensure(params && spectrum, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(cpl_table_get_nrow(spectrum) != 0, CPL_ERROR_DATA_NOT_FOUND, NULL);

    cpl_table *spec_telluriccorr_format = cpl_table_duplicate(spectrum);

    /* Initialize variables */
    cpl_error_code err = CPL_ERROR_NONE;

    /* Convert air to vacuum wavelengths if required */
    if (!err) {
        err = mf_spectrum_air_to_vacuum(params->config->inputs.wavelengths_frame, params->config->inputs.observing_erv_rv.value, spec_telluriccorr_format);
    }

    /* Save read parameter tables before manipulation */
    cpl_table *rangetab = NULL;
    cpl_table *chiptab  = NULL;
    if (!err) {
        rangetab = cpl_table_duplicate(params->rangetab);
        chiptab  = cpl_table_duplicate(params->chiptab );
    }

    cpl_size nchips = cpl_table_get_nrow(params->chiptab);



    /* Fill the parameter tables dependent on chip number */
    if (!err) err = mf_spectrum_fill_parameter_table(params, spec_telluriccorr_format);

    
    /* The following section is relevant to molecfit_model only. If this function has been */
    /* called by molecfit_calctrans then the supplied table inc)wrange will be NULL        */
    if (inc_wranges) {

	/* We are now at a point where we can check that range wavlengths ranges and */
	/* associated chip wavlength ranges are compatible                           */
	cpl_size nranges = cpl_table_get_nrow(inc_wranges);

	cpl_msg_info(cpl_func,"Nranges=%lld Nchips= %lld",nranges,nchips);

	/* Check if table inc_wranges has a chip mapping column */
	cpl_boolean has_chip_map = cpl_table_has_column(inc_wranges,MF_COL_WAVE_RANGE_MAP2CHIP);

	/* If no chip map then assert so and that all ranges will be mapped to chip 1 by default */
	if (!has_chip_map) {
            cpl_msg_info(cpl_func,"Include wavelength ranges table does not have a chip map.");
            cpl_msg_info(cpl_func,"Will assume that all ranges are mapped to chip 1");
	}

	/* Load up all the necessary chip and range wavelength associations for sanity checks */
	double chip_wave_beg[nchips];
	double chip_wave_end[nchips];
	double range_wave_beg[nranges];
	double range_wave_end[nranges];
	int    range_map2chip[nranges];
	
	
	for (cpl_size i=0; i<nchips; i++) {
	   chip_wave_beg[i]=cpl_table_get_double(params->chiptab,MF_COL_WL_MIN,i,NULL);
	   chip_wave_end[i]=cpl_table_get_double(params->chiptab,MF_COL_WL_MAX,i,NULL);
	}
	for (cpl_size i=0; i<nranges; i++) {
	   range_wave_beg[i]=cpl_table_get_double (inc_wranges,MF_COL_WAVE_RANGE_LOWER,   i,NULL);
	   range_wave_end[i]=cpl_table_get_double (inc_wranges,MF_COL_WAVE_RANGE_UPPER,   i,NULL);
	   if ( has_chip_map ) {
               range_map2chip[i]=cpl_table_get_int(inc_wranges,MF_COL_WAVE_RANGE_MAP2CHIP,i,NULL);
	   } else {
               range_map2chip[i]=1;
	   }
	}

   
	/* Count the number of chips specified in the map2chip vector */
	cpl_size n_mapped_chips = 0;
	for (cpl_size i=0; i<nranges; i++) {
            if (range_map2chip[i]>n_mapped_chips) n_mapped_chips=range_map2chip[i];
        }


        /* If there is only one chip being used for this run but there are many mapped chips then */
        /* not all ranges are being used and the following checks will be erroneous               */
        
        if (nchips==n_mapped_chips) {

	    /* Check that range wavelengths is contained within chip wavlength coverage */
	    for (cpl_size i=0; i<nranges; i++) {
	       int    chip_idx = range_map2chip[i]-1; 
	       double cwave_beg=chip_wave_beg[chip_idx];
	       double cwave_end=chip_wave_end[chip_idx];

	       /* Check the wavlength limits and ensure that there is a non zero overlap */
	       cpl_boolean OK_FLAG=CPL_TRUE;
	       if (range_wave_beg[i]>cwave_end) OK_FLAG=CPL_FALSE;
	       if (range_wave_end[i]<cwave_beg) OK_FLAG=CPL_FALSE;

	       /* Assert check results */
	       if (OK_FLAG) { 
                   cpl_msg_info(cpl_func,"Range %lld : [%f , %f] Mapped to Chip %d : [%f , %f] Overlap OK", 
	           i+1, range_wave_beg[i],range_wave_end[i],range_map2chip[i],cwave_beg,cwave_end);
	       } else {
                   cpl_msg_error(cpl_func,"Range %lld : %f , %f Mapped to Chip %d : %f , %f NO Overlaps. Not OK!",
	           i+1, range_wave_beg[i],range_wave_end[i],range_map2chip[i],cwave_beg,cwave_end);       
                   err=CPL_ERROR_ILLEGAL_INPUT;
	       } /* endif */

            } /* end for */
        
        } /* endif (nchips==n_mapped_chips) */
        

    } else {
    
        cpl_msg_info(cpl_func,"This is a single chip run amd some ranges are being mapped to chips not being used in this run");
    
    } /* end if (inc_wranges) */

    if (!err) {

        /* Default mrange column: full spectrum selected and range number for each chip */
        /*cpl_size nrow = cpl_table_get_nrow(inc_wranges);*/
        cpl_size nrow = cpl_table_get_nrow(spec_telluriccorr_format);
        for (cpl_size i = 0; i < nrow; i++) {
            int chip = cpl_table_get(spec_telluriccorr_format, MF_COL_CHIP,      i, NULL);
            cpl_table_set(           spec_telluriccorr_format, MF_COL_MOD_RANGE, i, chip);
        }

        /* Wavelength ranges for the fit */
        if (     !inc_wranges                        ) cpl_msg_info(cpl_func, "(mf_spectrum  ) Fit full spectrum (not provided inc_wranges fit ranges)");
        else if (cpl_table_get_nrow(inc_wranges) == 0) cpl_msg_info(cpl_func, "(mf_spectrum  ) Fit full spectrum (not provided inc_wranges fit ranges)");
        else {
            err = mf_spectrum_include_ranges(params, spec_telluriccorr_format, inc_wranges);
        }
    }


    /* Replace default parameter values in the range and chip tables by data from the parameter file */
    if (!err) err = mf_spectrum_replace_coefficients(params, rangetab, chiptab);

    if (!err && exc_wranges) err += mf_spectrum_exclude_ranges(spec_telluriccorr_format, exc_wranges, 'w');
    if (!err && exc_pranges) err += mf_spectrum_exclude_ranges(spec_telluriccorr_format, exc_pranges, 'p');

    /* Modify ranges to calculate only a single model spectrum? */
    if (!err) {

        /* Find largest step in wavenumber for range pixels */
        cpl_size nrow    = cpl_table_get_nrow(spec_telluriccorr_format);
        double   own     = 0.;
        double   maxstep = 0.;
        for (cpl_size i = 0; i < nrow; i++) {

            int range = cpl_table_get(spec_telluriccorr_format, MF_COL_MOD_RANGE, i, NULL);
            if (range > 0) {

                double nwn = 1e4 / cpl_table_get(spec_telluriccorr_format, MF_COL_IN_LAMBDA, i, NULL);
                double dwn = own - nwn;

                if (own > 0 && dwn > maxstep) maxstep = dwn;

                own = nwn;
            }
        }

        /* Criterion for single model spectrum */
        cpl_msg_info(cpl_func,"MNBX Criterion for single model spectrum %f, %d",maxstep,MF_EXTRA_WN_COVERAGE);
        if (maxstep > 2 * MF_EXTRA_WN_COVERAGE*0.1) {
            params->config->internal.single_spectrum = CPL_FALSE;
            cpl_msg_info(cpl_func,"MNBX FALSE");
        } else {
            cpl_msg_info(cpl_func,"MNBX TRUE");
            params->config->internal.single_spectrum = CPL_TRUE;

            /* Set all obsolete parameters to -1 */
            int nrange = params->config->internal.n_range;
            for (cpl_size i = 1; i < nrange; i++) {
                cpl_table_set(params->rangetab, MF_COL_PIX_RES,  i, -1);
                cpl_table_set(params->rangetab, MF_COL_WN_START, i, -1);
                cpl_table_set(params->rangetab, MF_COL_WN_END,   i, -1);
                cpl_table_set(params->rangetab, MF_COL_WN_STEP,  i, -1);
            }
        }
    }

    /* Set weights for spectra without error column */
    if (!err) err = mf_spectrum_set_weights(params, spec_telluriccorr_format);

    /* Check weights, 0? */
    if (!err && cpl_table_get_column_max(spec_telluriccorr_format, MF_COL_WEIGHT) == 0) err = CPL_ERROR_ILLEGAL_INPUT;

    /* Calculate required wavenumber ranges and store it in mf_parameters structure */
    if (!err) err = mf_spectrum_calculate_wave_range(params, spec_telluriccorr_format);

    if (!err && exc_wranges) err += mf_spectrum_exclude_ranges(spec_telluriccorr_format, exc_wranges, 'w');
    if (!err && exc_pranges) err += mf_spectrum_exclude_ranges(spec_telluriccorr_format, exc_pranges, 'p');


    /* Cleanup */
    if (rangetab) cpl_table_delete(rangetab);
    if (chiptab ) cpl_table_delete(chiptab );





    /* Check result and return */
    if (!err) {
        return spec_telluriccorr_format;
    } else {
        if (spec_telluriccorr_format) cpl_table_delete(spec_telluriccorr_format);
        return NULL;
    }
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
static cpl_error_code mf_spectrum_include_ranges(
    mf_parameters            *params,
    cpl_table                *spec,
    const cpl_table          *ranges)
{
    /*!
     * Reads wavelength ranges (in \f$\mu{\rm m}\f$) from a two-column ASCII
     * or FITS file (lower and upper limits) and selects these ranges in the
     * CPL table of the observed spectrum for the fitting procedure. Since the
     * radiative transfer code is only run for these ranges, the selection is
     * crucial for the code execution time. If a range table is not available,
     * the parameter 'wrange_include' in the mf_parameters structure has to be set
     * to 'none'. Apart from pairs of wavelength limits, an ASCII file can
     * also contain empty lines or comment lines starting with '#'.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param params  mf_parameters parameter structure
     *
     * \b OUTPUT:
     * \param spec    spectrum with mrange > 0 for selected ranges
     * \param params  mf_parameters with modified continuum fit parameters
     *
     * \b ERRORS:
     * - Unexpected file structure
     * - File opening failed
     * - see subroutines
     */

    /* Remove ranges outside valid wavelength range */
    double minlam = cpl_table_get_column_min(spec, MF_COL_IN_LAMBDA);
    double maxlam = cpl_table_get_column_max(spec, MF_COL_IN_LAMBDA);

    cpl_table *ranges_inc = cpl_table_duplicate(ranges);

    cpl_table_unselect_all(ranges_inc);
    cpl_table_or_selected_double(ranges_inc, MF_COL_WAVE_RANGE_LOWER, CPL_GREATER_THAN, maxlam);
    cpl_table_or_selected_double(ranges_inc, MF_COL_WAVE_RANGE_UPPER, CPL_LESS_THAN,    minlam);
    cpl_table_erase_selected(ranges_inc);

    /* Check number of remaining ranges */
    cpl_size nrange = cpl_table_get_nrow(ranges_inc);

    /* Create array for chip numbers of selected ranges */
    cpl_array *rangechips = cpl_array_new(nrange, CPL_TYPE_INT);

    /* Initialize column for range flag in spectrum table */
    cpl_size nrow = cpl_table_get_nrow(spec);
    cpl_table_fill_column_window_int(spec, MF_COL_MOD_RANGE, 0, nrow, 0);

    /* Set flags for ranges listed in the temporary CPL table */
    int range     = 0;
    int nrangeext = nrange;
    cpl_boolean isvalidranges = CPL_TRUE;
    for (cpl_size j = 0; j < nrange; j++) {

        double llim = cpl_table_get(ranges_inc, MF_COL_WAVE_RANGE_LOWER, j, NULL);
        double ulim = cpl_table_get(ranges_inc, MF_COL_WAVE_RANGE_UPPER, j, NULL);

        if (ulim <= llim) {
            isvalidranges = CPL_FALSE;
            nrangeext--;
            continue;
        }
        range++;

        cpl_size chip0 = -1;
        for (cpl_size i = 0; i < nrow; i++) {

            int    chip = cpl_table_get(spec, MF_COL_CHIP,      i, NULL);
            double lam  = cpl_table_get(spec, MF_COL_IN_LAMBDA, i, NULL);

            if (lam >= llim && lam <= ulim) {

                if (chip0 == -1) {
                    chip0 = cpl_table_get(spec, MF_COL_CHIP, i, NULL);
                }

                if (chip0 != chip) {
                    /* New range number if chip changes */
                    range++;
                    nrangeext++;
                    chip0 = chip;
                    cpl_array_set_size(rangechips, nrangeext);
                }

                cpl_table_set(spec, MF_COL_MOD_RANGE, i, range);
                cpl_array_set(rangechips, range-1, chip0);
            }
        }
    }

    /* Set number of ranges in parameter list */
    params->config->internal.n_range = nrangeext;

    /* Get fit flag for continuum fit */
    cpl_boolean fit_cont = params->config->fitting.fit_continuum.fit;

    /* Expand table of continuum fit parameters in mf_parameters structure */
    cpl_table_set_size(params->rangetab, nrangeext);
    cpl_array *cont_coef = cpl_array_duplicate(cpl_table_get_array(params->rangetab, MF_COL_CONT_COEF, 0));

    for (cpl_size j = 0; j < nrangeext; j++) {
        cpl_table_set(      params->rangetab, MF_COL_CHIP,      j, cpl_array_get(rangechips, j, NULL));
        cpl_table_set(      params->rangetab, MF_COL_FIT_RANGE, j, fit_cont                          );
        cpl_table_set_array(params->rangetab, MF_COL_CONT_COEF, j, cont_coef                         );
    }

    /* Cleanup */
    cpl_table_delete(ranges_inc);
    cpl_array_delete(rangechips);
    cpl_array_delete(cont_coef);

    /* Management of errors */
    if (isvalidranges == CPL_FALSE) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                                     "Unexpected cpl_table structure: Upper limit <= lower limit");
    }

    cpl_table_and_selected_double(spec, MF_COL_WEIGHT, CPL_GREATER_THAN, 0.);
    int n_weight_selected = cpl_table_count_selected(spec);
    cpl_table_select_all(spec);
    if (n_weight_selected == 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): cpl_table *spec (all weights = 0)");
    }

    cpl_table_and_selected_int(spec, MF_COL_MOD_RANGE, CPL_GREATER_THAN, 0);
    int n_ranges_selected = cpl_table_count_selected(spec);
    cpl_table_select_all(spec);
    if (nrange > 0 && n_ranges_selected == 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): cpl_table *spec (wavelength inclusion range selects no values in the input data)");
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
static cpl_error_code mf_spectrum_exclude_ranges(
    cpl_table                *spec,
    const cpl_table          *ranges,
    const char               wp)
{
    /*!
     * Reads wavelength (in \f$\mu{\rm m}\f$) or pixel ranges (indicated by
     * 'w' or 'p' for input parameter 'wp') from a two-column ASCII or FITS
     * file (lower and upper limits) and sets the weight of these ranges in
     * the CPL table of the observed spectrum to zero. If a range table is not
     * available, the parameter 'wrange_exclude' or 'prange_exclude' in the
     * mf_parameters structure has to be set to 'none'. Apart from pairs of
     * wavelength or pixel limits, an ASCII file can also contain empty lines
     * or comment lines starting with '#'.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param params  mf_parameters parameter structure
     * \param wp      w[avelength] or p[ixel] ranges
     *
     * \b OUTPUT:
     * \param spec    spectrum with weight = 0 for excluded ranges
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     * - Unexpected file structure
     * - File opening failed
     * - see subroutines
     */

    /* Get file type of range table */
    cpl_msg_info(cpl_func, "(mf_spectrum  ) mf_spectrum_exclude_ranges(...)");

    cpl_size    nrow          = cpl_table_get_nrow(spec);
    double      val           = 0.;
    cpl_boolean isvalidranges = CPL_TRUE;
    int         nrange        = cpl_table_get_nrow(ranges);

    /* Set weight to zero for ranges listed in the temporary CPL table */
    for (cpl_size j = 0; j < nrange && isvalidranges; j++) {

        double llim = cpl_table_get(ranges, MF_COL_WAVE_RANGE_LOWER, j, NULL);
        double ulim = cpl_table_get(ranges, MF_COL_WAVE_RANGE_UPPER, j, NULL);

        if (ulim < llim) {

            isvalidranges = CPL_FALSE;

        } else {

            for (cpl_size i = 0; i < nrow; i++) {

                if (     wp == 'w') val = cpl_table_get(spec, MF_COL_IN_LAMBDA, i, NULL);   /* Wavelength */
                else if (wp == 'p') val = i + 1;                                            /* Pixel      */

                if (val >= llim && val <= ulim) cpl_table_set(spec, MF_COL_WEIGHT, i, 0.);
            }
        }
    }

    /* Error management */
    if (isvalidranges == CPL_FALSE) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                                     "Unexpected cpl_table structure: Upper limit < lower limit");
    }

    cpl_table_and_selected_double(spec, MF_COL_WEIGHT, CPL_GREATER_THAN, 0.);
    double n_weight_selected = cpl_table_count_selected(spec);
    cpl_table_select_all(spec);
    if (n_weight_selected == 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): cpl_table *spec (exclusion range rejects all values in the input data)");
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
static cpl_error_code mf_spectrum_replace_coefficients(
    mf_parameters            *params,
    cpl_table                *rangetab,
    cpl_table                *chiptab)
{
    /*!
     * Replaces entries in the range and chip tables of the mf_parameters parameter
     * structure by values provided by the parameter file. In the range table,
     * the range-specific fit flags and the coefficients of the continuum fit
     * are updated. In the chip table, the chip-specific fit flags and the
     * coefficients of the wavelength correction are modified. Only those
     * table entries are changed which were explicitly given in the parameter
     * file using the keywords fit_range[range], cont_range[range],
     * fit_chip[chip], and wlc_chip[chip].
     *
     * \b INPUT:
     * \param params    mf_parameters parameter structure
     * \param rangetab  range table with data from the input parameter file
     * \param chiptab   chip table with data from the input parameter file
     *
     * \b OUTPUT:
     * \param params    mf_parameters with modified range and chip tables
     *
     * \b ERRORS:
     * - none
     */

   
    /* Get the number of physical wavelength ranges */
    cpl_size nrange = cpl_table_get_nrow(params->rangetab);
      

    /* Adapt size of input range table with data to be added */
    cpl_table_set_size(rangetab, nrange);

    cpl_msg_info(cpl_func,"Number of physical wavelength ranges=%lld",nrange);

    /* Adapt maximum number of continuum coefficients in input range table with data to be added */
    cpl_size ncont = cpl_table_get_column_depth(params->rangetab, MF_COL_CONT_COEF);
    cpl_table_set_column_depth(rangetab, MF_COL_CONT_COEF, ncont);

    /* Replace range-specific flags for continuum fit if they were set in the parameter file */
    for (cpl_size i = 0; i < nrange; i++) {
        if (cpl_table_is_valid( rangetab, MF_COL_FIT_RANGE, i) == 1) {
            int fit = cpl_table_get(rangetab, MF_COL_FIT_RANGE, i, NULL);
            cpl_table_set(params->rangetab,   MF_COL_FIT_RANGE, i, fit );
	    cpl_msg_info(cpl_func,"STORED FLAG FOR RANGE %lld is valid and set to %d (not substituting)",i+1,fit);
        } else {
	    int fit=0;
	    if (params->config->fitting.fit_ranges[i]) fit=1;
	    cpl_msg_info(cpl_func,"No external definition of fit flag for range %lld Parameter flag set to %d",i+1,fit);
            cpl_table_set(params->rangetab, MF_COL_FIT_RANGE, i, fit);
	}
    }



    /* ------------------------------------------------------------------------------*/
    /* Initial values for the coefficients of the RANGE-SPECIFIC continuum fit models*/
    /* ------------------------------------------------------------------------------*/
    if (params->config->fitting.expert_mode) {
        cpl_msg_info(cpl_func,"EXPERT MODE FLAG=TRUE so using coefficient values for continuum model that were imported");
        cpl_msg_info(cpl_func,"No of ranges = %lld", nrange);
    
        /* We populate the coefficients array elements of the chip table from values
	   stored in params->config->fitting.cnt_coeffs[][]   */
           
           cpl_table_dump(params->rangetab,0,nrange,NULL);
    
	for (int i=0;i<nrange;i++) {

            int cnt_porder=params->config->fitting.cont_poly_order[i];
            int max_porder=params->config->fitting.fit_continuum.n;
            cpl_msg_info(cpl_func,"Poly order for the continuum model of range %d  is %d",i+1,cnt_porder);
	    cpl_array *coef_out = cpl_array_new(max_porder+1,CPL_TYPE_DOUBLE);
            for (cpl_size j = 0; j <=max_porder; j++) cpl_array_set(coef_out, j, 0.0);

            for (cpl_size j = 0; j <=cnt_porder; j++) {
                    double coef =params->config->fitting.cont_coeffs[i][j];
                    cpl_array_set(coef_out, j, coef);
		    cpl_msg_info(cpl_func,"Replace CONTINUUM coefs for range %d coef %lld new coef=%f",i+1,j,coef);
            }/*end j*/
	
            cpl_table_set_array(params->rangetab, MF_COL_CONT_COEF, i, coef_out);
	    cpl_array_delete(coef_out);
	
	}/* end i*/

    } else {

        cpl_msg_info(cpl_func,"EXPERT MODE FLAG=FALSE");

        /* Replace range-specific continuum fit coefficients if they were set in the parameter file */
        for (cpl_size i = 0; i < nrange; i++) {

            if (!cpl_table_is_valid(rangetab, MF_COL_CONT_COEF, i)) {
                continue;
            }

            cpl_array *coef_in  = cpl_array_duplicate(cpl_table_get_array(rangetab,         MF_COL_CONT_COEF, i));
            cpl_array *coef_out = cpl_array_duplicate(cpl_table_get_array(params->rangetab, MF_COL_CONT_COEF, i));

            for (cpl_size j = 0; j < ncont; j++) {
                if (cpl_array_is_valid(coef_in, j) == 1) {
                    double coef = cpl_array_get(coef_in, j, NULL);
                    cpl_array_set(coef_out, j, coef);
                }
            }

            cpl_table_set_array(params->rangetab, MF_COL_CONT_COEF, i, coef_out);

            cpl_array_delete(coef_in);
            cpl_array_delete(coef_out);
        }
    
    } /*endif*/

    /* Adapt size of input chip table with data to be added */
    cpl_size nchip = cpl_table_get_nrow(params->chiptab);
    cpl_table_set_size(chiptab, nchip);
 
    /* Adapt maximum number of wavelength coefficients in input chip table with data to be added */
    cpl_size nwlc = cpl_table_get_column_depth(params->chiptab, MF_COL_WLC_COEF);
    cpl_table_set_column_depth(chiptab, MF_COL_WLC_COEF, nwlc);

    cpl_msg_info(cpl_func,"nchip=%lld. Polynomial order=%lld",nchip,nwlc-1);


    /* Replace chip-specific flags for wavelength fit if they were set in the parameter file */
    for (cpl_size i = 0; i < nchip; i++) {
        if (cpl_table_is_valid( chiptab, MF_COL_FIT_CHIP, i) == 1) {
            int fit = cpl_table_get(chiptab, MF_COL_FIT_CHIP, i, NULL);
            cpl_table_set(params->chiptab, MF_COL_FIT_CHIP, i, fit);
        } else {
	    int fit=0;
	    if (params->config->fitting.fit_chips[i]) fit=1;
	    cpl_msg_info(cpl_func,"No external definition of fit flag for chip %lld Parameter flag set to %d",i+1,fit);
            cpl_table_set(params->chiptab, MF_COL_FIT_CHIP, i, fit);
 	}
    }


    /* ------------------------------------------------------------------------------*/
    /* Initial values for the coefficients of the CHIP-SPECIFIC wavelength fit models*/
    /* ------------------------------------------------------------------------------*/
    if (params->config->fitting.expert_mode) {
        cpl_msg_info(cpl_func,"EXPERT MODE FLAG=TRUE");
    
        /* We populate the coefficients array elements of the chip table from values
	   stored in params->config->fitting.wlc_coeffs[][]   */
    
	int    wlc_porder      = params->config->fitting.fit_wavelenght.n;
	for (int i=0;i<nchip;i++) {

            cpl_array *coef_out = cpl_array_new(wlc_porder+1,CPL_TYPE_DOUBLE);

            for (cpl_size j = 0; j <=wlc_porder; j++) {
                    double coef =params->config->fitting.wlc_coeffs[i][j];
                    cpl_array_set(coef_out, j, coef);
		    cpl_msg_info(cpl_func,"Replace WLC coefs for chip %d coef %lld new coef=%f",i+1,j,coef);
            }/*end j*/
	
            cpl_table_set_array(params->chiptab, MF_COL_WLC_COEF, i, coef_out);
	    cpl_array_delete(coef_out);
	
	}/* end i*/

    } else {

        cpl_msg_info(cpl_func,"EXPERT MODE FLAG=FALSE");
        /* The original code pre "EXPERT_MODE" implementation */

	for (cpl_size i = 0; i < cpl_table_get_nrow(chiptab); i++) {

            if (!cpl_table_is_valid(chiptab, MF_COL_WLC_COEF, i)) {
        	continue;
            }

            cpl_array *coef_in  = cpl_array_duplicate(cpl_table_get_array(chiptab,         MF_COL_WLC_COEF, i));
            cpl_array *coef_out = cpl_array_duplicate(cpl_table_get_array(params->chiptab, MF_COL_WLC_COEF, i));

            for (cpl_size j = 0; j < nwlc; j++) {
        	if (cpl_array_is_valid(coef_in, j) == 1) {
                    double coef = cpl_array_get(coef_in, j, NULL);
                    cpl_array_set(coef_out, j, coef);
       	        }
            }

            cpl_table_set_array(params->chiptab, MF_COL_WLC_COEF, i, coef_out);

            cpl_array_delete(coef_in);
            cpl_array_delete(coef_out);
	}

    }/* endf expert_mode*/

    
    return CPL_ERROR_NONE;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief Convert wavelengths to vacuum.
 *
 * @param frame_type         Type of wavelength reference.
 * @param obs_erf_rv         The observatory radial velocity relative o an external reference frame (typically, the sun or the barycenter of the solar system)
 * @param spectrum           cpl_table with the spectrum (in microns).
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
static cpl_error_code mf_spectrum_air_to_vacuum(
    const char               *frame_type,
    const double             obs_erf_rv,
    cpl_table                *spectrum)
{
    /* No action is performed if the wavelengths of the input spectrum are already for vacuum.*/
    if (!strcmp(frame_type, MF_PARAMETERS_WAVELENGTH_FRAME_VACUUM)) return CPL_ERROR_NONE;

    /* Get number of rows and wavelengths */
    cpl_size nrow    = cpl_table_get_nrow(spectrum);
    double   *lambda = cpl_table_get_data_double(spectrum, MF_COL_IN_LAMBDA);

    if (!strcmp(frame_type, MF_PARAMETERS_WAVELENGTH_FRAME_VACUUM_RV)) {

        /* The input wavelength has been corrected for the Earth motion relative to an external reference frame */
        for (cpl_size i = 0; i < nrow; i++) {

            double n = (1. + 1.55e-8) * (1. + obs_erf_rv / 299792.458);

            lambda[i] /= n;
        }

    } else if (!strcmp(frame_type, MF_PARAMETERS_WAVELENGTH_FRAME_AIR)) {

        /* Converts air wavelengths to vacuum wavelengths by using the formula of Edlen (1966) : Get refractive index and convert wavelengths. */
        for (cpl_size i = 0; i < nrow; i++) {

             double sig2 = pow(lambda[i], -2);
             double n    = 8342.13 + 2406030. / (130. - sig2) + 15997. / (38.9 - sig2);

             n = 1. + 1e-8 * n;

             lambda[i] *= n;
        }

    } else {

        /* Illegal input frame_type value */
        return CPL_ERROR_ILLEGAL_INPUT;
    }

    return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Fills the parameter tables of the mf_parameters parameter structure for the number of chips.
 *
 * @param params    Parameters structure.
 * @param spectrum  cpl_table with observed spectrum.
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK. cpl_error_code in other cases
 *
 * @description This function adapt the parameters structure to the number of chips
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_spectrum_fill_parameter_table(
    mf_parameters            *params,
    cpl_table                *spectrum)
{
   /*!
     * Fills the parameter tables of the mf_parameters parameter structure for the number of chips.
     *
     * \b INPUT:
     * \param params  mf_parameters parameter structure
     * \param spec    CPL table with observed spectrum
     *
     * \b OUTPUT:
     * \param params  mf_parameters with parameter tables adapted to number of chips
     *
     * \b ERRORS:
     * - none
     */

    /* Get number of chips and polynomial coefficients */
    cpl_size nwlc  = 1 + params->config->fitting.fit_wavelenght.n;
    cpl_size nchip = params->config->internal.nchip;
    cpl_msg_info(cpl_func,"fit_wavelength n=%lld, nwlc=%lld", params->config->fitting.fit_wavelenght.n,nwlc);
    
    /* Note that the parameter params->config->fitting.fit_continuum.n is the maximum */
    /*  poly order of all the ranges                                                  */
    cpl_size ncont = 1 + params->config->fitting.fit_continuum.n;

    

    /* Declare only 1 range per chip */
    params->config->internal.n_range = nchip;

    /* Get fit flag for continuum fit */
    cpl_boolean fit_cont = params->config->fitting.fit_continuum.fit;

    /* Resize table for range-related parameters */
    cpl_table_set_size(params->rangetab, nchip);
    cpl_table_set_column_depth(params->rangetab, MF_COL_CONT_COEF, ncont);


    /*** Fill table for range-related parameters ***/

    /* Coefficients for continuum correction */
    cpl_array *cont_coef = cpl_array_new(ncont, CPL_TYPE_DOUBLE);
    cpl_array_fill_window(cont_coef, 0, ncont, 0.);

    cpl_array_set(cont_coef, 0, params->config->fitting.fit_continuum.const_val);
    for (cpl_size chip_index = 0; chip_index < nchip; chip_index++) {

        cpl_size chip = chip_index + 1;

        cpl_table_set(      params->rangetab, MF_COL_CHIP,      chip_index, chip     );
        cpl_table_set(      params->rangetab, MF_COL_FIT_RANGE, chip_index, fit_cont );
        cpl_table_set_array(params->rangetab, MF_COL_CONT_COEF, chip_index, cont_coef);
    }
    cpl_array_delete(cont_coef);

    /* Get fit flag for wavelength correction */
    cpl_boolean fit_wlc = params->config->fitting.fit_wavelenght.fit;

    /* Get constant term for wavelength correction */
    double wlc0 = params->config->fitting.fit_wavelenght.const_val;

    /* Resize table for chip-related parameters */
    cpl_table_set_size(params->chiptab, nchip);
    cpl_table_set_column_depth(params->chiptab, MF_COL_WLC_COEF, nwlc);


    /*** Fill table for chip-related parameters ***/

    /* Coefficients for wavelength correction */
    cpl_array *wlc_coef = cpl_array_new(nwlc, CPL_TYPE_DOUBLE);
    cpl_array_fill_window(wlc_coef, 0, nwlc, 0.);
    cpl_array_set(wlc_coef, 0, wlc0);
    cpl_array_set(wlc_coef, 1, 1.);
    cpl_table_fill_column_window(      params->chiptab, MF_COL_FIT_CHIP, 0, nchip, fit_wlc );
    cpl_table_fill_column_window_array(params->chiptab, MF_COL_WLC_COEF, 0, nchip, wlc_coef);
    cpl_array_delete(wlc_coef);

    /* Minimum and maximum wavelengths */

    for (cpl_size chip_index = 0; chip_index < nchip; chip_index++) {

        cpl_size chip = chip_index + 1;

        /* Extract wavelength grid of observed spectrum for each chip */
        cpl_table_unselect_all(spectrum);
        cpl_table_or_selected_int(spectrum, MF_COL_CHIP, CPL_EQUAL_TO, chip);
        cpl_table *chip_spec_selected = cpl_table_extract_selected(spectrum);
        cpl_table_select_all(spectrum);

        cpl_size nsel = cpl_table_get_nrow(chip_spec_selected);
        cpl_msg_info(cpl_func,"--> NSEL =%lld CHIP_INDEX=%lld",nsel,chip_index);

        /* Get minimum and maximum wavelengths and write them into table */
        cpl_table_set(params->chiptab, MF_COL_WL_MIN, chip_index, cpl_table_get(chip_spec_selected, MF_COL_IN_LAMBDA, 0,        NULL));
        cpl_table_set(params->chiptab, MF_COL_WL_MAX, chip_index, cpl_table_get(chip_spec_selected, MF_COL_IN_LAMBDA, nsel - 1, NULL));

        /* Delete temporary table */
        cpl_table_delete(chip_spec_selected);
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
static cpl_error_code mf_spectrum_set_weights(
    const mf_parameters      *params,
    cpl_table                *spec)
{
    /*!
     * Sets the weights for spectra that do not contain an error column.
     * The default relative error provided by the mf_parameters parameter structure
     * is multiplied by the mean flux of all wavelengths that are used by the
     * fit. The resulting absolute error is taken for all wavelengths with
     * non-zero weight.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param params  mf_parameters parameter structure
     *
     * \b OUTPUT:
     * \param spec    spectrum with optimised weights
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    /* Return if error data exist */
    if (strcmp(params->config->inputs.column_dflux, MF_PARAMETERS_COLUMN_DFLUX_NULL) != 0) {
        return CPL_ERROR_NONE;
    }

    /* Get default relative error */
    double deferr = params->config->inputs.default_error;

    /* Count data points with non-zero weight and sum up their fluxes */
    int      nw   = 0;
    double   fsum = 0.;
    cpl_size nrow = cpl_table_get_nrow(spec);
    for (cpl_size i = 0; i < nrow; i++) {

        double weight = cpl_table_get(spec, MF_COL_WEIGHT, i, NULL);

        if (weight > 0) {
            nw++;
            fsum += cpl_table_get(spec, MF_COL_IN_FLUX, i, NULL);
        }
    }

    /* Set errors to default error * mean if error column is missing */
    if (fsum > 0) {

        double weight0 = (double)nw / (fsum * deferr);

        for (cpl_size i = 0; i < nrow; i++) {
            double weight = cpl_table_get(spec, MF_COL_WEIGHT, i, NULL);
            if (weight > 0) cpl_table_set(spec, MF_COL_WEIGHT, i, weight0);
        }
    }

    /* Check for errors */
    cpl_error_code err = CPL_ERROR_NONE;
    if (     fsum == 0) err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                                    "Invalid object value(s): cpl_table *spec (all fluxes = 0)");
    else if (nw   == 0) err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                                    "Invalid object value(s): cpl_table *spec (all weights = 0)");
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
static cpl_error_code mf_spectrum_calculate_wave_range(
    mf_parameters            *params,
    cpl_table                *spec)
{
    /*!
     * Calculates the model wavenumber limits and step sizes for the different
     * selected wavelength ranges and stores this information in the mf_parameters
     * parameter structure.
     *
     * \b INPUT:
     * \param params  mf_parameters parameter structure
     * \param spec    CPL table with observed spectrum
     *
     * \b OUTPUT:
     * \param params  mf_parameters parameter structure with wavenumber limits and
     *                and step sizes
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - No data
     */

    /* Get slit width in pixels -> rough estimate/lower limit of resolution in pixels */
    double slit_width  = params->config->instrumental.slit_width.value;
    double pixel_scale = params->config->instrumental.pixel_scale.value;
    if (pixel_scale == 0.) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid object value(s): pixel_scale of mf_parameters *params = 0");
    }

    /* Require at least 2 pixels per resolution element */
    double respix = slit_width / pixel_scale;
    if (respix < 2.) respix = 2.;

    /* Check existence of data */
    int m = cpl_table_get_nrow(spec);
    if (m <= 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "No data: cpl_table *spec");
    }

    /* Get number of wavelength ranges : Calculate one molecular spectrum for full wavelength range? */
    int nrange = params->config->internal.n_range;
    if (params->config->internal.single_spectrum) {
        nrange = 1;
    }

    /* Calculate input parameters for radiative transfer code */
    for (int j = 0; j < nrange; j++) {

        /* Extract wavelength grid of observed spectrum for each range */
        cpl_table_unselect_all(spec);
        if (params->config->internal.single_spectrum) {
            cpl_table_or_selected_int(spec, MF_COL_MOD_RANGE, CPL_NOT_EQUAL_TO, 0    );
        } else {
            cpl_table_or_selected_int(spec, MF_COL_MOD_RANGE, CPL_EQUAL_TO,     j + 1);
        }

        cpl_table *rangespec = cpl_table_extract_selected(spec);
        cpl_table_select_all(spec);

        cpl_size nsel = cpl_table_get_nrow(rangespec);

        /* Case: all weights = 0 */
        if (cpl_table_get_column_max(rangespec, MF_COL_WEIGHT) == 0.) {
            cpl_msg_info(cpl_func, "(mf_spectrum  ) Range %d: all weights = 0", j + 1);
            cpl_table_set(params->rangetab, MF_COL_PIX_RES,  j, 0.);
            cpl_table_set(params->rangetab, MF_COL_WN_START, j, 0.);
            cpl_table_set(params->rangetab, MF_COL_WN_END,   j, 0.);
            cpl_table_set(params->rangetab, MF_COL_WN_STEP,  j, 0.);
            cpl_table_delete(rangespec);
            continue;
        }

        /* Derive wavenumber range in cm^-1 */
        double wn_start = 1e4 / cpl_table_get(rangespec, MF_COL_IN_LAMBDA, nsel - 1, NULL);
        double wn_end   = 1e4 / cpl_table_get(rangespec, MF_COL_IN_LAMBDA, 0,        NULL);

        /* Consider additional spectral coverage */
        wn_start -= MF_EXTRA_WN_COVERAGE;
        wn_end   += MF_EXTRA_WN_COVERAGE;

        /* Find highest resolution */
        double hres = 0.;
        double olam = 0.;
        for (cpl_size i = 0; i < nsel; i++) {
            double nlam = cpl_table_get(rangespec, MF_COL_IN_LAMBDA, i, NULL);
            double res  = nlam / ((nlam - olam) * respix);
            if (res > hres) {
              hres = res;
            }
            olam = nlam;
        }

        /* Derive wavenumber step */
        double dwn     = wn_end - wn_start;
        double wn_step = wn_end / (hres * respix * MF_OVERSAMPLING_FACTOR);
        double npix    = dwn / wn_step;

        wn_step *= npix / ceil(npix);

        /* Write pixel resolution, wavenumber range and step into mf_parameters structure */
        cpl_table_set(params->rangetab, MF_COL_PIX_RES,  j, hres * respix);
        cpl_table_set(params->rangetab, MF_COL_WN_START, j, wn_start     );
        cpl_table_set(params->rangetab, MF_COL_WN_END,   j, wn_end       );
        cpl_table_set(params->rangetab, MF_COL_WN_STEP,  j, wn_step      );

        /* Delete temporary table */
        cpl_table_delete(rangespec);
    }

   return CPL_ERROR_NONE;
}

/** @endcond */


/**@}*/
