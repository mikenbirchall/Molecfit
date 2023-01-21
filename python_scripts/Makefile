ODA_DIR=OpticalDepthAdaptation
BIN_DIR=bin

all: $(BIN_DIR) ODA_scripts

ODA_scripts: lblrtmOda lnflOda mbmOda mrgOda

lblrtmOda:
	cp $(ODA_DIR)/$@.py $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

lnflOda:
	cp $(ODA_DIR)/$@.py $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

mbmOda:
	cp $(ODA_DIR)/$@.py $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

mrgOda:
	cp $(ODA_DIR)/$@.py $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

$(BIN_DIR):
	mkdir $@

clean:
	rm -rdf $(BIN_DIR)
