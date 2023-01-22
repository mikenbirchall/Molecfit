ODA_DIR=OpticalDepthAdaptation
BIN_DIR=bin

all: $(BIN_DIR) ODA_scripts

ODA_scripts: lblrtmOda lnflOda mbmOda mrgOda plot28

lblrtmOda:
	cat $(ODA_DIR)/CommonUtils.py $(ODA_DIR)/$@.py > $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

lnflOda:
	cat $(ODA_DIR)/CommonUtils.py $(ODA_DIR)/$@.py > $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

mbmOda:
	cat $(ODA_DIR)/CommonUtils.py $(ODA_DIR)/$@.py > $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

mrgOda:
	cat $(ODA_DIR)/CommonUtils.py $(ODA_DIR)/$@.py > $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

plot28:
	cat $(ODA_DIR)/CommonUtils.py $(ODA_DIR)/$@.py > $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

$(BIN_DIR):
	mkdir $@

clean:
	rm -rdf $(BIN_DIR)
