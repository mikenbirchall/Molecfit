HTTPS_PATH = https://ftp.emc.ncep.noaa.gov/jcsda/CRTM/utility/lblrtmio/v0.1.0
top_dir= $(shell pwd)
lblrtmio_ver = lblrtmio_v0.1.0
lblrtmio_dir = $(top_dir)/$(lblrtmio_ver)
tarball      = $(lblrtmio_ver).tar.gz
curlpath     = $(HTTPS_PATH)/$(tarball)
build_dir    = $(top_dir)/build
bin_dir      = $(build_dir)/bin
lib_dir      = $(build_dir)/$(lblrtmio_ver)/lib
lblrtmio_arv = $(lib_dir)/liblblrtmio.a
src_dir      = $(top_dir)/src
CAT_FILE     = $(src_dir)/cat.txt
EXE_FILES =$(shell cat $(CAT_FILE) | tr '\n' ' ')

all: $(build_dir) lblrtmio binaries


lblrtmio: $(build_dir) $(lblrtmio_dir)
	cd $(lblrtmio_dir); ./configure --prefix=$(build_dir)\
	                                --enable-single
	make all install -C $(lblrtmio_dir)

$(build_dir):
	if [ ! -e $(build_dir) ] ; then  mkdir $(build_dir) ; fi
	if [ ! -e $(bin_dir)   ] ; then  mkdir $(bin_dir)    ; fi

$(lblrtmio_dir): $(tarball)
	tar -zxvf $(tarball)

$(tarball):
	curl $(HTTPS_PATH)/$(tarball) -o $(tarball)

clean:
	cd src; make clean

cleandist: clean
	rm -rdf $(lblrtmio_dir)
	rm -rdf $(build_dir)

binaries:
	make -C $(src_dir) all
	echo $(EXE_FILES)

install:
	if [ -e $(CAT_FILE) ] ; then cd $(src_dir) ; cp $(EXE_FILES) $(bin_dir) ; fi
