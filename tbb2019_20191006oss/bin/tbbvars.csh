#!/bin/csh
#
# Copyright (c) 2005-2019 Intel Corporation
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set LIBTBB_NAME="libtbb.dylib"

# Arg1 represents TBBROOT detection method. Its possible value is 'auto_tbbroot'. In which case
# the environment variable TBBROOT is detected automatically by using the script directory path.
if ("$1" == "auto_tbbroot") then
    set sourced=($_)
    if ("$sourced" != '') then # if the script was sourced
        set tmp="`dirname "$sourced[2]"`"
        set script_dir=`cd "$tmp"; pwd -P`
    else # if the script was run => "$_" is empty
        echo "ERROR: you have to source the script."
        exit 1
    endif
    setenv TBBROOT "$script_dir/.."
else
    setenv TBBROOT "SUBSTITUTE_INSTALL_DIR_HERE"
endif

if (-e "${TBBROOT}/lib/${LIBTBB_NAME}") then
    if (! $?LIBRARY_PATH) then 
        setenv LIBRARY_PATH "${TBBROOT}/lib"
    else
        setenv LIBRARY_PATH "${TBBROOT}/lib:$LIBRARY_PATH"
    endif
    if (! $?DYLD_LIBRARY_PATH) then
        setenv DYLD_LIBRARY_PATH "${TBBROOT}/lib"
    else
        setenv DYLD_LIBRARY_PATH "${TBBROOT}/lib:$DYLD_LIBRARY_PATH"
    endif
    if (! $?CPATH) then
        setenv CPATH "${TBBROOT}/include"
    else
        setenv CPATH "${TBBROOT}/include:$CPATH"
    endif
else
    echo "ERROR: ${LIBTBB_NAME} library does not exist in ${TBBROOT}/lib."
    unsetenv TBBROOT
    exit 1
endif
