#!/bin/bash
cd ${0%/*} || exit 1    # run from this directory

export CURRDIR="$PWD"
export SRC=$CURRDIR/src
export FO=$SRC/functionObjects

# Function Objects
wmake libso $FO/wettedArea
wmake libso $FO/contactAngleEvaluation

# Utilities

#######################################



