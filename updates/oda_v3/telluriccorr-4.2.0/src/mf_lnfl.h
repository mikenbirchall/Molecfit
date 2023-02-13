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

#ifndef MF_LNFL_H
#define MF_LNFL_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>

#include "mf_constants.h"
#include "mf_io.h"
#include "mf_parameters.h"

CPL_BEGIN_DECLS

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

/*** line database for the preparation with LNFL ***/

/* The aer line list is delivered with the LNFL/LBLRTM package [1]. In the case of an update, the new file name has to be adapted.
 * In addition, there might be a change in the file format. Please see [1] for more details.
 *
 * Caveat: The aer line list is optimised for the usage with LNFL/LBLRTM.
 *         A change of the incorporated file or the format might lead to unexcpected and unreliable results.
 *
 *         [1] http://rtweb.aer.com/lblrtm_frame.html
 */
#define MF_LNFL_LINE_DB                  "LNFL_LINE_DB"
#define MF_LNFL_LINE_DB_DESC             "File name of the line list (must be stored in the directory :\n"\
                                         " ({TELLURICCORR_DATA_PATH}/hitran/)"
#define MF_LNFL_LINE_DB_AER_VERSION      "aer_v_3.8.1.2"
#define MF_LNFL_LINE_DB_INIT             MF_LNFL_LINE_DB_AER_VERSION

/* Format of the line file: gives the length in terms of characters per line
 *
 * Is usually = 100 (old HITRAN format) or = 160 (new HITRAN format) */
#define MF_LNFL_LINE_DB_FMT              "LNFL_LINE_DB_FORMAT"
#define MF_LNFL_LINE_DB_FMT_DESC         "Format of the line file: gives the length in terms of characters per line"
#define MF_LNFL_LINE_DB_FMT_OLD_HITRAN   100
#define MF_LNFL_LINE_DB_FMT_NEW_HITRAN   160
#define MF_LNFL_LINE_DB_FMT_INIT         MF_LNFL_LINE_DB_FMT_OLD_HITRAN

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

/* Create the user parameter structure (with default values) */
MF_EXPORT mf_io_lnfl_config * mf_lnfl_config_create(void);

/* Deallocate the user parameter structure */
MF_EXPORT void mf_lnfl_config_delete(
    mf_io_lnfl_config        *config);        /* In/out: parameter structure  */

/* Check the user parameter structure. For check the user modifications */
MF_EXPORT cpl_error_code mf_lnfl_config_check(
    mf_io_lnfl_config        *config);        /* User parameter structure     */

/*  */
MF_INTERNAL cpl_error_code mf_lnfl(
    mf_io_lnfl_config        *config,
    mf_parameters            *params);


CPL_END_DECLS


#endif /* MF_LNFL_H */
