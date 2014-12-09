cmake_minimum_required (VERSION 2.6)
project ( lime )
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG -Wall ")

find_package(PythonLibs 2.7 REQUIRED)

file(GLOB limeSRC
    "./include/*.h"
    "./src/*.c"
)

include_directories(include ${PYTHON_INCLUDE_DIRS} )
add_executable( lime ${limeSRC} )
target_link_libraries( lime -lgsl -lm -lcfitsio -lcblas -lcurses -lqhull ${PYTHON_LIBRARIES} )


