# Please ensure that any changes remain compliant with 2.4.8.
if(NOT EMBED_OPENBABEL)
  cmake_minimum_required(VERSION 2.4.8)
endif()
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )
SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/test" )
set(OPENBABEL_BINARY_ROOT "/home/cyx/workspace/ob/")
set(PROJECT_SOURCE_DIR "/home/cyx/workspace/DEPACT_release")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -D_GLIBCXX_DEBUG -fno-inline -Wl,--enable-new-dtags -Wl,--fatal-warnings -Wl,--no-undefined -lc -ldl -lz")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -L/usr/local/lib/ -lmlpack -lblas -lboost_serialization -lboost_filesystem -lboost_system")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O4 -std=c++11 -Wl,--no-as-needed -L${PROJECT_BINARY_DIR}/lib/libcombine.a -L/usr/local/lib/ -lboost_serialization -lboost_filesystem -lboost_system -lgslcblas -lgsl -L${PROJECT_SOURCE_DIR}/cliquer/libcliquer.a")
# Quick: -O4
link_directories(${PROJECT_BINARY_DIR}/lib  ${OPENBABEL_BINARY_ROOT}/lib  ${PROJECT_SOURCE_DIR}/cliquer)
include_directories(${PROJECT_SOURCE_DIR} ${OPENBABEL_BINARY_ROOT}/include/openbabel-2.0  ${PROJECT_SOURCE_DIR}/noob)
add_library(combine STATIC IMPORTED GLOBAL)
set(C_LIB ${CMAKE_BINARY_DIR}/lib/libcombine.a)
set_target_properties(combine
	PROPERTIES
	IMPORTED_LOCATION ${C_LIB}
	)
add_dependencies(combine combinelib)
add_custom_target(combinelib
	COMMAND rm -f libcombine.a
	COMMAND ar rcT libcombine.a libmycode.a libnoob.a libnspdataio.a libnspgeometry.a  libnspproteinrep.a
        COMMAND ranlib libcombine.a
	DEPENDS mycode noob nspdataio nspproteinrep cliquer
        WORKING_DIRECTORY ${LIBRARY_OUTPUT_PATH})

add_subdirectory(src)
add_subdirectory(noob)
add_subdirectory(test)
add_subdirectory (dataio/src)
add_subdirectory (dstl/src)
add_subdirectory (geometry/src)
add_subdirectory (proteinrep/src)
