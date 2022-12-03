PWD_DIR = $(shell pwd)
TOP_DIR = $(shell dirname $(PWD_DIR))
LIB_DIR = $(TOP_DIR)/build/lblrtmio_v0.1.0/lib
INC_DIR = $(TOP_DIR)/build/lblrtmio_v0.1.0/include
ARV_LIB = $(LIB_DIR)/liblblrtmio.a
CAT_FILE  = cat.txt
EXE_FILES =$(shell cat $(CAT_FILE) | tr '\n' ' ')
TARGETS = test ODdefltParse

list :
	@echo $(EXE_FILES)

all: clean $(TARGETS)

ODdefltParse:
	make  generic TARGET=$@

test:
	gfortran -c $@.f90 -J$(INC_DIR)
	gfortran  $@.o -J$(INC_DIR) $(ARV_LIB) -o $@
	@echo Built executable $@
	echo $@ >> $(CAT_FILE)

generic:
	gfortran -c $(TARGET).f90 -J$(INC_DIR)
	gfortran  $(TARGET).o -J$(INC_DIR) $(ARV_LIB) -o $(TARGET)
	@echo Built executable $(TARGET)
	echo $(TARGET) >> $(CAT_FILE)

clean:
	rm -f *.o
	rm -f $(TARGETS)
	rm -f $(CAT_FILE)

.PHONY: generic
