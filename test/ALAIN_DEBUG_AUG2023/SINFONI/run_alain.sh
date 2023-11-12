#!/usr/bin/tcsh

if ($?MOLECFIT_DATA_HOME) then
    echo Enviromental variable MOLECFIT_DATA_HOME is set to $MOLECFIT_DATA_HOME
else
    echo Enviromental variable MOLECFIT_DATA_HOME is not defined
    exit
endif
if ($?MOLECFIT_INSTALL_DIR) then
    echo Enviromental variable MOLECFIT_INSTALL_DIR is set to $MOLECFIT_INSTALL_DIR
else
    echo Enviromental variable MOLECFIT_INSTALL_HOME is not defined
    exit
endif

setenv BIN_DIR ${MOLECFIT_INSTALL_DIR}/bin
setenv PATH    ${BIN_DIR}:${PATH}

setenv SOF_DATA ${MOLECFIT_DATA_HOME}/ALAIN_DEBUG_AUG2023/SINFONI
rm -rdf output; mkdir output
rm -rdf logs  ; mkdir logs
esorex --recipe-config=config/model.rc --output-dir=output --log-dir=logs molecfit_model sof/model.sof
