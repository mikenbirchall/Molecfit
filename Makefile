token =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
github=https:///github.com/mikenbirchall/Molecfit

molecfit-kit-ver=molecfit-kit-4.2.3
tarball = $(molecfit-kit-ver).tar.gz
https://ftp.eso.org/pub/dfs/pipelines/instruments/molecfit/molecfit-kit-4.2.3.tar.gz
tarball_url=https://ftp.eso.org/pub/dfs/pipelines/instruments/molecfit

TOP_DIR=molecfit-kit-4.2.3
PWD=$(shell pwd)
IDR=$(PWD)/install
LIB_DIR=$(IDR)/lib

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
	@echo "    all          ---> Builds for standard kit"
	@echo "    updates      ---> Updates standrad kit source code from specifc update"
	@echo "    rebuild      ---> Rebuilds, eg after updated changes"
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
	@echo "    scripts      ---> Creates shell scipts to setup PATH variable to execute binaries"

all: decompress cfitsio fftw wcslib cpl esorex third_party telluriccorr molecfit scripts

tarball:
	make $(tarball)

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

telluriccorr:
	cd $(telluriccorr_dir); ./configure --prefix=$(IDR) --with-cpl=$(IDR)
	cd $(telluriccorr_dir); make all install

molecfit:
	cd $(molecfit_dir); ./configure --prefix=$(IDR) --with-cpl=$(IDR)\
							--with-telluriccorr=$(IDR)
	cd $(molecfit_dir); make all install

scripts:
	@echo export PATH=$(IDR):$$PATH > setup.sh
	@echo setenv PATH $(IDR):$$PATH > setup.csh
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
	rm -rdf install
	mkdir install

cleankit:
	rm -rdf $(molecfit-kit-ver)


updates:
	@echo Updating $(TOP_DIR) kit from: updates/$(UPDATES)
	@make update4target UPDATE_TARGET=$(telluriccorr_dir)/src

update4target:
# This recipe does the bulk of the file copy routine
# for updating the kit files and relies on the same subdirectory
# structure being used in the target and source directores.
	@echo Updating source files in directory: $(UPDATE_TARGET)
	@cp ${subst $(TOP_DIR),$(UPDATES_DIR),$(UPDATE_TARGET)/*} \
	$(UPDATE_TARGET)

rebuild:
	cd $(telluriccorr_dir); make all install

token:
	@cat $(HOME)/github.dat


.PHONY: updates
