
# This file is part of the LITIV framework; visit the original repository at
# https://github.com/plstcharles/litiv for more information.
#
# Copyright 2018 Pierre-Luc St-Charles; pierre-luc.st-charles<at>polymtl.ca
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

if(CMAKE_VERSION VERSION_LESS 3.2)
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
else()
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
endif()

if((NOT DEFINED GEOS_SOURCE_DIR) OR (NOT EXISTS ${GEOS_SOURCE_DIR}))
    download_project(
        PROJ
            GEOS
        PREFIX
            "${CMAKE_CURRENT_BINARY_DIR}"
        URL
            "http://download.osgeo.org/geos/geos-3.8.1.tar.bz2"
        QUIET
    #    ${UPDATE_DISCONNECTED_IF_AVAILABLE} # bugged in some cmake versions, needs patch
    )
    set(GEOS_SOURCE_DIR ${GEOS_SOURCE_DIR} CACHE INTERNAL "")
    set(GEOS_BINARY_DIR ${GEOS_BINARY_DIR} CACHE INTERNAL "")
endif()

list(
	APPEND
		CMAKE_MODULE_PATH
	"${GEOS_SOURCE_DIR}/cmake/"
)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

add_subdirectory(
	${GEOS_SOURCE_DIR}
	${GEOS_BINARY_DIR}
	EXCLUDE_FROM_ALL
)
