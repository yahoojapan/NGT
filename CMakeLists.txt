cmake_minimum_required(VERSION 2.8)

project( ngt )

set(ngt_VERSION_MAJOR 0 )
set(ngt_VERSION_MINOR 4 )
set(ngt_VERSION_PATCH 0 )

set( ngt_VERSION ${ngt_VERSION_MAJOR}.${ngt_VERSION_MINOR}.${ngt_VERSION_PATCH} )
set( ngt_SOVERSION ${ngt_VERSION_MAJOR} )

if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if( ${UNIX} )
	include_directories("/home/y/include64/" "/home/y/include/")
	link_directories("/home/y/lib64" "/usr/lib64")

	set(CMAKE_SKIP_BUILD_RPATH TRUE)
	set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
	set(CMAKE_INSTALL_RPATH "/home/y/lib64")

	set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -std=gnu++0x") # for std::unordered_set
	set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O2 -std=gnu++0x") # for std::unordered_set

	add_subdirectory("${PROJECT_SOURCE_DIR}/lib")
	add_subdirectory("${PROJECT_SOURCE_DIR}/bin")
endif( ${UNIX} )
