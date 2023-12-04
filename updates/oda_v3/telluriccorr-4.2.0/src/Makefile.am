## Process this file with automake to produce Makefile.in

##   This file is part of the ESO Telluric Correction Library
##   Copyright (C) 2001-2018 European Southern Observatory

##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.

##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111-1307  USA

AUTOMAKE_OPTIONS = 1.8 foreign

SUBDIRS = . test

AM_CFLAGS = \
            -Werror \
            $(OPENMP_CFLAGS)

AM_CPPFLAGS = \
              -DTELLURICCORR_BASEDIR='"$(prefix)"' \
              -DTELLURICCORR_DATADIR='"$(profilesdir)"' \
              $(CPL_INCLUDES)

DISTCLEANFILES = \
                 *~ \
                 *.o \
                 *.out \
                 *.gcno \
                 *.gcda \
                 Makefile.in

noinst_includedir =
noinst_include_HEADERS = \
                         cpl_mpfit.h \
                         mf_constants.h \
                         mf_io.h \
                         mf_strlst.h \
                         mf_kernel_user.h \
                         mf_kernel_synthetic.h \
                         mf_parameters.h \
                         mf_spectrum.h \
                         mf_gdas.h \
                         mf_atm_combined_gdas.h \
                         mf_fit.h \
                         mf_molecules.h \
                         mf_lnfl.h \
                         mf_lblrtm.h \
                         mf_convolution.h \
                         mf_configuration.h \
                         mf_model.h \
                         mf_calctrans.h \
                         mf_wrap.h \
                         mf_wrap_utils.h \
                         mf_wrap_data.h \
                         mf_wrap_config.h \
                         mf_wrap_fits.h \
                         mf_wrap_dfs.h \
                         mf_wrap_calc.h \
                         mf_wrap_config_calc.h \
                         mf_wrap_config_corr.h

                         

include_HEADERS = \
                  ../telluriccorr.h \
                  mf_constants.h \
                  mf_io.h \
                  mf_strlst.h \
                  mf_parameters.h \
                  mf_kernel_user.h \
                  mf_molecules.h \
                  mf_spectrum.h \
                  mf_gdas.h \
                  mf_atm_combined_gdas.h \
                  mf_lnfl.h \
                  mf_lblrtm.h \
                  mf_convolution.h \
                  mf_configuration.h \
                  mf_model.h \
                  mf_calctrans.h \
                  mf_wrap.h \
                  mf_wrap_utils.h \
                  mf_wrap_data.h \
                  mf_wrap_config.h \
                  mf_wrap_fits.h \
                  mf_wrap_dfs.h \
                  mf_wrap_calc.h \
                  mf_wrap_config_calc.h \
                  mf_wrap_config_corr.h \
                  defs.h \
                  cpp.h \
                  c.h \
                  cerf.h


lib_LTLIBRARIES = libtelluriccorr.la

libtelluriccorr_la_SOURCES = \
                         cpl_mpfit.c \
                         mf_io.c \
                         mf_strlst.c \
                         mf_parameters.c \
                         mf_kernel_user.c \
                         mf_kernel_synthetic.c \
                         mf_spectrum.c \
                         mf_gdas.c \
                         mf_atm_combined_gdas.c \
                         mf_fit.c \
                         mf_molecules.c \
                         mf_lnfl.c \
                         mf_lblrtm.c \
                         mf_convolution.c \
                         mf_configuration.c \
                         mf_model.c \
                         mf_calctrans.c \
                         mf_wrap.c \
                         mf_wrap_utils.c \
                         mf_wrap_data.c \
                         mf_wrap_config.c \
                         mf_wrap_fits.c \
                         mf_wrap_dfs.c \
                         mf_wrap_calc.c  \
                         mf_wrap_config_calc.c \
                         mf_wrap_config_corr.c \
                         w_of_z.c \
                         width.c \
                         im_w_of_x.c\
                         err_fcts.c\
                         erfcx.c
libtelluriccorr_la_CFLAGS = $(AM_CFLAGS)
libtelluriccorr_la_LDFLAGS = $(CPL_LDFLAGS) -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libtelluriccorr_la_LIBADD = $(LIBCPLDFS) $(LIBCPLUI) $(LIBCPLDRS) $(LIBCPLCORE)
