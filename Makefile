github=https:///github.com/mikenbirchall/Molecfit

molecfit-kit-ver=molecfit-kit-4.2.3
tarball = $(molecfit-kit-ver).tar.gz
tarball_url=https://ftp.eso.org/pub/dfs/pipelines/instruments/molecfit

# GSL STUFF
gsl_url    =https://ftp.gnu.org/gnu/gsl
gsl_ver    =gsl-2.5
gsl_tarball=$(gsl_ver).tar.gz
gsl_home   =gsl
gsl_dir    =$(gsl_home)/$(gsl_ver)


TOP_DIR=molecfit-kit-4.2.3
PWD=$(shell pwd)
IDR=$(PWD)/install
LIB_DIR=$(IDR)/lib
BIN_DIR=$(IDR)/bin
CPLDR=$(IDR)

cfitsio_dir     =$(TOP_DIR)/cfitsio-3.49
wcslib_dir      =$(TOP_DIR)/wcslib-7.6
fftw_dir        =$(TOP_DIR)/fftw-3.3.9
cpl_dir         =$(TOP_DIR)/cpl-7.2.2
esorex_dir      =$(TOP_DIR)/esorex-3.13.6
thirdparty_dir  =$(TOP_DIR)/molecfit_third_party-1.9.2
telluriccorr_dir=$(TOP_DIR)/telluriccorr-4.2.0
molecfit_dir    =$(TOP_DIR)/molecfit-4.2.3

UPDATES=continuum_parameters
UPDATES_DIR=updates/$(UPDATES)

what:
	@echo "    all          ---> Builds Molecfit based on standard kit and updated with default ODA version"
	@echo "    noesorex     ---> As make all, but doesnt build esorex. Must specify CPLDR eg make esorex CPLDR=[path ro cpl]"
	@echo "    esorex       ---> Just builds esorex from standard kit"
	@echo "    base         ---> Builds base Molecfit from standard kit (no upddates added)"
	@echo "    updates      ---> Updates standard kit source code from specifc update, eg, make updates UPDATES=continuum_parameters"
	@echo "    rebuild      ---> Rebuilds, eg after updated changes"
	@echo "    rebuildwgsl  ---> Rebuilds with the wgsl mod and the gsl library. Note gsl library must be built first"
	@echo "    allwgsl      ---> Special case of Build for standard kit with gsl library and updates from UPDATES=wgsl"
	@echo "    tarball      ---> Fetch, if need be, the molecfit kit tarball"
	@echo "    decompress   ---> Decompres the molecfit kit tarball"
	@echo "    cfitsio      ---> Build the cfitsio component"
	@echo "    fftw         ---> Build the fftw component"
	@echo "    wcslib       ---> Build the wcslib component"
	@echo "    cpl          ---> Build the cpl component"
	@echo "    esorex       ---> Build the esorex component"
	@echo "    third_party  ---> Build the third_party component"
	@echo "    telluriccorr ---> Build the telluricor component"
	@echo "    molecfit     ---> Build the molecfit recipe component"
	@echo "    lblrtmio     ---> Builds lblrtmio binaries"
	@echo "    python_scripts--> Install python scripts"
	@echo "    scripts      ---> Creates shell scripts to setup PATH variable to execute binaries"
	@echo "    gslib        ---> Adds, builds and installs the GNU Scientific library to the install directory"
	@echo "    lmods        ---> List all update mods available to build with"


all: alloda

base: decompress cfitsio fftw wcslib cpl esorex third_party telluriccorr molecfit

kitver:
	make base IDR=${PWD}/install-${molecfit-kit-ver}

allwgsl: base
	make gslib
	make updates UPDATES=wgsl
	make rebuildwgsl

alloda: base oda

esorexonly: decompress cfitsio fftw wcslib cpl esorex

# If not building esorex but pointing to another installation with CPLDR i.e.
#         make noesorex CPLDR=[path to CPL installation]
# This recipe checks the validity of CPLDR
chk_cpldr:
	@echo Checking specified CPL directory...
	@if [   -d $(CPLDR) ];  then echo Specified CPL Directory Found ; fi
	@if [ ! -d $(CPLDR) ];  then echo Specified CPL Directory $(CPLDR) not found. Aborting! && exit 1 ; fi
	@echo Checking that CPL directory is different from the installation directory...
ifeq ($(IDR) , $(CPLDR))
	@echo CPL Directory specified is the same as the installation directory, Error assumed.
	@echo Will abort
	@exit 1
else
	@echo Directories are different will continue
endif
	@echo Check that specified CPL directory contains the cpl core lib...
	@if [ -e $(CPLDR)/lib/libcplcore.la ]; then echo cplcore found ; else echo cplcore not found. Aborting && exit 1 ; fi

# Quick route to build oda from scratch without building esorex but specifying a
# location for a prebuilt CPL directory (CPLDR)
noesorex: chk_cpldr decompress third_party telluriccorr molecfit oda

oda:
	make updates UPDATES=oda_v3
	make rebuild
	make python_scripts

tarball:
	# If no tarball then make (fetch) one
	if [ ! -e $(tarball) ] ; then make $(tarball); fi


$(tarball):
	curl $(tarball_url)/$(tarball) -o $(tarball)

decompress: cleankit $(tarball)
	tar -zxvf $(molecfit-kit-ver).tar.gz
	cd $(TOP_DIR); tar -zxvf cfitsio-3.49.tar.gz
	cd $(TOP_DIR); tar -zxvf fftw-3.3.9.tar.gz
	cd $(TOP_DIR); bunzip2   wcslib-7.6.tar.bz2
	cd $(TOP_DIR); tar -xvf  wcslib-7.6.tar
	cd $(TOP_DIR); tar -zxvf cpl-7.2.2.tar.gz
	cd $(TOP_DIR); tar -zxvf esorex-3.13.6.tar.gz
	cd $(TOP_DIR); tar -xvf  molecfit_third_party-1.9.2.tar
	cd $(TOP_DIR); tar -zxvf telluriccorr-4.2.0.tar.gz
	cd $(TOP_DIR); tar -zxvf molecfit-4.2.3.tar.gz

cfitsio:
	cd $(cfitsio_dir); ./configure --prefix=$(IDR) --enable-shared --enable-reentrant
	cd $(cfitsio_dir); make all
	cd $(cfitsio_dir); make shared
	cd $(cfitsio_dir); make install

wcslib:
	cd $(wcslib_dir); ./configure --prefix=$(IDR) --with-cfitsio=$(IDR)
	cd $(wcslib_dir); make all install

fftw:
	cd $(fftw_dir); ./configure --prefix=$(IDR) --enable-shared   \
	                                            --disable-fortran \
	                                            --enable-sse2
	cd $(fftw_dir); make all install clean
	cd $(fftw_dir); ./configure --prefix=$(IDR) --enable-float    \
	                                            --enable-shared   \
	                                            --disable-fortran \
	                                            --enable-sse      \
	                                            --enable-sse2
	cd $(fftw_dir); make all install
cpl:
	cd $(cpl_dir); export LD_LIBRARY_PATH=$(LIB_DIR);\
		       ./configure --prefix=$(IDR) --with-cfitsio=$(IDR) \
	                                           --with-wcslib=$(IDR)  \
	                                           --with-fftw=$(IDR)
	cd $(cpl_dir); make all install

esorex:
	cd $(esorex_dir); ./configure --prefix=$(IDR) --with-cpl=$(IDR) \
	                                           --with-wcslib=$(IDR)
	cd $(esorex_dir); make all install

third_party:
	cd $(thirdparty_dir); make -f BuildThirdParty.mk all install prefix=$(IDR)

third_party_dbl:
	cd $(thirdparty_dir); make -f BuildThirdParty.mk clean all \
				      lnfl_target=customLinuxGNUdbl
	cp $(thirdparty_dir)/lblrtm/lblrtm*dbl $(BIN_DIR)/lblrtm_dbl


telluriccorr:
	cd $(telluriccorr_dir); ./configure --prefix=$(IDR) --with-cpl=$(CPLDR)
	cd $(telluriccorr_dir); make all install

molecfit:
	cd $(molecfit_dir); ./configure --prefix=$(IDR) --with-cpl=$(CPLDR)\
							--with-telluriccorr=$(IDR)
	cd $(molecfit_dir); make all install


gslib:
	# If no gsl home directory then make one
	if [ ! -e $(gsl_home) ] ; then mkdir $(gsl_home); fi

	# If no gsl tarball in home directory then curl one
	if [ ! -e $(gsl_home)/$(gsl_tarball) ] ; then \
	    cd $(gsl_home); curl $(gsl_url)/$(gsl_tarball) -o $(gsl_tarball); \
	fi

	# If no untarred toplevel gsl directory then untar one
	if [ ! -e $(gsl_dir) ] ; then \
	    cd $(gsl_home); tar -zxvf $(gsl_tarball); \
	fi

	# If no makefile in the gsl top directory then use configure
	# with specification of where the installation directory is.
	if [ ! -e $(gsl_dir)/Makefile ] ; then \
	    cd $(gsl_dir); ./configure --prefix=$(IDR); \
	fi

	# Now make all and install
	cd $(gsl_dir); make all install

scripts:
	@echo export PATH=$(IDR)/bin:$$PATH > setup.sh
	@echo setenv PATH $(IDR)/bin:$$PATH > setup.csh
	@echo Shell scripts to include built bin to PATH var created.
	@echo From bash shell type:   "source setup.sh"
	@echo From c/tcsh shell type: "source setup.csh"

clean_scripts:
	rm -f setup.sh
	rm -f setup.csh

distclean: clean_scripts cleaninstall
	cd $(cfitsio_dir); if [ -e Makefile ] ; then make distclean ; fi
	cd $(wcslib_dir);  if [ -e Makefile ] ; then make distclean ; fi
	cd $(fftw_dir);    if [ -e Makefile ] ; then make distclean ; fi
	cd $(cpl_dir);     if [ -e Makefile ] ; then make distclean ; fi

cleaninstall:
	rm -rdf $(IDR)
	mkdir $(IDR)

cleankit:
	rm -rdf $(molecfit-kit-ver)

cleanall: distclean cleankit

updates:
	@echo Updating $(TOP_DIR) kit from: updates/$(UPDATES)
	@make update4target UPDATE_TARGET=$(telluriccorr_dir)/src

update4target:
# This recipe does the bulk of the file copy routine
# for updating the kit files and relies on the same subdirectory
# structure being used in the target and source directores.
	@echo Updating source files in directory: $(UPDATE_TARGET)
	@cp ${subst $(TOP_DIR),$(UPDATES_DIR),$(UPDATE_TARGET)/*.c} $(UPDATE_TARGET)
	@cp ${subst $(TOP_DIR),$(UPDATES_DIR),$(UPDATE_TARGET)/*.h} $(UPDATE_TARGET)

rebuild:
	cd $(telluriccorr_dir); make all install

rebuildwgsl:
	# First Copy the wgsl mods
	make updates UPDATES=wgsl
	#To rebuild with gsl we have to clean out the src directory then call make with specific flag variable
	cd $(telluriccorr_dir)/src; make clean
	cd $(telluriccorr_dir)/src; make all libtelluriccorr_la_LIBADD="\$$(LIBCPLDFS) \$$(LIBCPLUI) \$$(LIBCPLDRS) \$$(LIBCPLCORE) -lgsl -lgslcblas -lm"
	cd $(telluriccorr_dir)/src; make install

lblrtmio:
	cd lblrtmio; make all install
	cp lblrtmio/build/bin/* $(BIN_DIR)

python_scripts:
	cd python_scripts; make all
	cp python_scripts/bin/* $(BIN_DIR)
	cd python_scripts; make clean

token:
	@cat $(HOME)/github.dat

lmods:
	cd perl_scripts; perl diff_search.pl ..

.PHONY: updates lblrtmio python_scripts
