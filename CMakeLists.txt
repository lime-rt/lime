cmake_minimum_required (VERSION 2.6)
project ( lime )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG")

file(GLOB limeSRC
    "./include/*.h"
    "./src/*.c"
)

include_directories(include)
add_executable( lime ${limeSRC} )
target_link_libraries( lime -lgsl -lm -lcfitsio -lcblas -lcurses -lqhull -lmatheval )


