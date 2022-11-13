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

#ifndef MF_CONVOLUTION_H
#define MF_CONVOLUTION_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/
#include "../../install/include/gsl/gsl_linalg.h"
#include "../../install/include/gsl/gsl_matrix.h"

#include <cpl.h>

#include "mf_constants.h"
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

#define MF_CONVOLUTION_LAM_UNIT       1e-6          /* mu m                   */
#define MF_CONVOLUTION_SR_IN_ARCSEC2  4.254517e+10  /* steradians in arcsec^2 */

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
MF_INTERNAL cpl_error_code mf_convolution(
    const mf_parameters      *params,
    const cpl_boolean        correct_spectrum,
    const cpl_boolean        last_call,
    cpl_table                **spec_out,
    cpl_error_code           *range_status,
    cpl_array                *fitpar,
    int                      mpfit_calls,
    cpl_table                *spec,
    cpl_matrix               *kernel);

CPL_END_DECLS


#endif /* MF_CONVOLUTION_H */
