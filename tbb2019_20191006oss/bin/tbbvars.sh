#!/bin/sh
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

LIBTBB_NAME="libtbb.dylib"
TBBROOT="SUBSTITUTE_INSTALL_DIR_HERE"
# Arg1 represents TBBROOT detection method. Its possible value is 'auto_tbbroot'. In which case
# the environment variable TBBROOT is detected automatically by using the script directory path.
if [ -n "${BASH_SOURCE}" ]; then
    if [ "$1" = "auto_tbbroot" ]; then
        TBBROOT=$(cd $(dirname ${BASH_SOURCE}) && pwd -P)/..
    fi
fi

if [ -e "$TBBROOT/lib/$LIBTBB_NAME" ]; then
    export TBBROOT
    if [ -z "${LIBRARY_PATH}" ]; then
        LIBRARY_PATH="${TBBROOT}/lib"; export LIBRARY_PATH
    else
        LIBRARY_PATH="${TBBROOT}/lib:$LIBRARY_PATH"; export LIBRARY_PATH
    fi
    if [ -z "${DYLD_LIBRARY_PATH}" ]; then
        DYLD_LIBRARY_PATH="${TBBROOT}/lib"; export DYLD_LIBRARY_PATH
    else
        DYLD_LIBRARY_PATH="${TBBROOT}/lib:$DYLD_LIBRARY_PATH"; export DYLD_LIBRARY_PATH
    fi
    if [ -z "${CPATH}" ]; then
        CPATH="${TBBROOT}/include"; export CPATH
    else
        CPATH="${TBBROOT}/include:$CPATH"; export CPATH
    fi
else
    echo "ERROR: $LIBTBB_NAME library does not exist in $TBBROOT/lib."
    return 1
fi
