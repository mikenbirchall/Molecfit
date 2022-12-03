PWD_DIR = $(shell pwd)
TOP_DIR = $(shell dirname $(PWD_DIR))
LIB_DIR = $(TOP_DIR)/build/lblrtmio_v0.1.0/lib
INC_DIR = $(TOP_DIR)/build/lblrtmio_v0.1.0/include
ARV_LIB = $(LIB_DIR)/liblblrtmio.a
CAT_FILE  = cat.txt
EXE_FILES =$(shell cat $(CAT_FILE) | tr '\n' ' ')
TARGETS = test

list :
	@echo $(EXE_FILES)

all: clean $(TARGETS)

test:
	gfortran -c $@.f90 -J$(INC_DIR)
	gfortran  $@.o -J$(INC_DIR) $(ARV_LIB) -o $@
	@echo Built executable $@
	echo $@ >> $(CAT_FILE)
	cp $@ test23
	echo test23 >> $(CAT_FILE)
	cp $@ test42
	echo test42 >> $(CAT_FILE)

list1 :
	echo $(PWD_DIR) $(TOP_DIR)
	ls $(LIB_DIR)
	ls $(INC_DIR)

clean:
	rm -f *.o
	rm -f $(TARGETS)
	rm -f $(CAT_FILE)
