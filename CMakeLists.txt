cmake_minimum_required (VERSION 2.6)
project ( lime )
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG -Wall ")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(PythonLibs 2.7 REQUIRED)
find_package(GSL REQUIRED)
find_package(M REQUIRED)
find_package(CFITSIO REQUIRED)
find_package(CBLAS REQUIRED)
find_package(Curses REQUIRED)
find_package(Qhull REQUIRED)

file(GLOB limeSRC
    "${CMAKE_SOURCE_DIR}/include/*.h"
    "${CMAKE_SOURCE_DIR}/src/*.c"
)

set(LIBS ${LIBS} ${PYTHON_LIBRARIES} ${GSL_LIBRARIES} ${M_LIBRARIES} ${CFITSIO_LIBRARIES} ${CBLAS_LIBRARIES} ${CURSES_LIBRARIES} ${QHULL_LIBRARIES} )
include_directories(include ${PYTHON_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${M_INCLUDE_DIRS} ${CFITSIO_INCLUDE_DIRS} ${CBLAS_INCLUDE_DIRS} ${CURSES_INCLUDE_DIRS} ${QHULL_INCLUDE_DIRS} )
configure_file( ./doc/example/input.ini.in ./doc/example/input.ini @ONLY )
add_executable( lime ${limeSRC} )
target_link_libraries( lime ${LIBS} )


