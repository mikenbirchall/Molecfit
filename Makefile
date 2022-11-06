molecfit-kit-ver=molecfit-kit-4.2.3

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
UPDATES_DIR=$(UPDATES)/$(telluriccorr_dir)


all: decompress cfitsio fftw wcslib cpl esorex third_party telluriccorr molecfit scripts

decompress: cleankit
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
	echo Updating from: $(UPDATES)
	ls $@/$(UPDATES)/$(telluriccorr_dir)/src
	cp $@/$(UPDATES)/$(telluriccorr_dir)/src/*.c $(telluriccorr_dir)/src
	cp $@/$(UPDATES)/$(telluriccorr_dir)/src/*.h $(telluriccorr_dir)/src
	cd $(telluriccorr_dir); make all install


.PHONY: updates
