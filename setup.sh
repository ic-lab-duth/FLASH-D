#! /bin/bash

#================================================================
# File: setup.sh
#
# Usage: source setup.sh
#
# This script will setup any dependency need by 
# this project.

# It will create local directories:
#    ./ac_types
#    ./ac_simutils
#    ./ac_math
#    ./prompt_bench
#================================================================

# Configure AC Datatypes
if [ ! -d ./ac_types ]; then
  echo "Downloading AC_Types..."
  git clone https://github.com/hlslibs/ac_types.git
fi
export AC_TYPES=$(pwd)/ac_types

# Configure AC Simutils
if [ ! -d ./ac_simutils ]; then
  echo "Downloading AC_Simutils..."
  git clone https://github.com/hlslibs/ac_simutils.git
fi
export AC_SIMUTILS=$(pwd)/ac_simutils

# Configure AC Math
if [ ! -d ./ac_math ]; then
  echo "Downloading AC_Math..."
  git clone https://github.com/hlslibs/ac_math.git
fi
export AC_SIMUTILS=$(pwd)/ac_math

# Configure Prompt Bench
if [ ! -d ./prompt_bench ]; then
  echo "Downloading PromptBench..."
  git clone https://github.com/microsoft/promptbench.git
fi
export PROMPT_BENCH=$(pwd)/prompt_bench

# TODO: Add python dependencies