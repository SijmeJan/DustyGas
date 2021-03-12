#!/bin/bash

if [ "$(uname)" == "Darwin" ]; then
    # Running on MacOS

    # Check if greadlink is installed
    command -v greadlink >/dev/null 2>&1 || { echo >&2 "Requires greadlink but it's not installed.  Aborting."; exit 1; }

    echo "Replacing readlink with greadlink in updateSubmodules"
    awk '{gsub("\\(readlink", "\(greadlink")}1' ExaHyPE-Engine/Submodules/updateSubmodules.sh > test.dat_tmp && chmod +x test.dat_tmp && mv test.dat_tmp ExaHyPE-Engine/Submodules/updateSubmodules.sh

    echo "Replacing readlink with greadlink in updateSubmodules"
    awk '{gsub("\\(readlink", "\(greadlink")}1' ExaHyPE-Engine/Toolkit/toolkit.sh > test.dat_tmp  && chmod +x test.dat_tmp && mv test.dat_tmp ExaHyPE-Engine/Toolkit/toolkit.sh

    echo "Setting SHAREDMEM=None"
    export SHAREDMEM=None

fi

ExaHyPE-Engine/Submodules/updateSubmodules.sh

if [ "$(uname)" == "Darwin" ]; then
    # Running on MacOS

    # Comment out all -lrt lines in Makefile
    awk '/lrt/{
        if (substr($1,1,1) !~ /^[#]/)
        {
                a = $0;
                sub(/[^\t].*/, "", a);
                $1 = "#" $1;
                $0 = a $0
        }
    }{print}' ExaHyPE-Engine/ExaHyPE/Makefile > makefile_tmp && mv makefile_tmp ExaHyPE-Engine/ExaHyPE/Makefile
fi

# Works on Linux, not on my Macbook
# jobs.h needs extra argument 'priority'?
# Core.cpp uses header file not available on MacOS: set cores to 1?
