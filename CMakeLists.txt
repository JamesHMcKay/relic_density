cmake_minimum_required(VERSION 2.8)
project(main)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(src include)
set(HEADER_FILES include/rd.hpp include/interp.hpp  include/integrate.hpp include/models.hpp include/bessel.hpp include/data.hpp)
set(SRC_FILES src/main.cpp src/singletdm.cpp src/singletdmz3.cpp src/interp.cpp src/rd.cpp )


set(Boost_INCLUDE_DIR /usr/local/Cellar/boost/1.59.0/include)
set(Boost_LIBRARY_DIR /usr/local/Cellar/boost/1.59.0/lib)
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})



add_executable(main ${SRC_FILES} ${HEADER_FILES})

target_link_libraries( main ${Boost_LIBRARIES} )


if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS main DESTINATION ${CMAKE_BINARY_DIR}/lib)