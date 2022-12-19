ODA_DIR=OpticalDepthAdaptation
BIN_DIR=bin

all: $(BIN_DIR) ODA_scripts

ODA_scripts: lblrtmOda lnflOda

lblrtmOda:
	cp $(ODA_DIR)/$@.py $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

lnflOda:
	cp $(ODA_DIR)/$@.py $(BIN_DIR)/$@
	chmod a+x $(BIN_DIR)/$@

$(BIN_DIR):
	mkdir $@

clean:
	rm -rdf $(BIN_DIR)
