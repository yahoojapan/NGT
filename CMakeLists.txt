cmake_minimum_required(VERSION 2.8)

project( ngt )

set(ngt_VERSION_MAJOR 1 )
set(ngt_VERSION_MINOR 2 )
set(ngt_VERSION_PATCH 1 )

set( ngt_VERSION ${ngt_VERSION_MAJOR}.${ngt_VERSION_MINOR}.${ngt_VERSION_PATCH} )
set( ngt_SOVERSION ${ngt_VERSION_MAJOR} )

if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if( ${UNIX} )
    if( CMAKE_VERSION VERSION_LESS 3.1 )
	link_directories("/usr/lib64")

	set(CMAKE_SKIP_BUILD_RPATH TRUE)
	set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

	if( ${NGT_AVX_DISABLED} )
	    	message(STATUS "AVX will not be used to compute distances.")
		set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -std=gnu++0x -lrt") # for std::unordered_set, std::unique_ptr
		set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3 -std=gnu++0x -lrt") # for std::unordered_set, std::unique_ptr
	else()
		set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -std=gnu++0x -lrt -mavx") # for std::unordered_set, std::unique_ptr
		set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3 -std=gnu++0x -lrt -mavx") # for std::unordered_set, std::unique_ptr
	endif()
    else()
  	option(WALL "enable all warnings" ON)
 	if( ${WALL} )
 		add_compile_options(-Wall)
 	endif()
 
 	# CMAKE_CXX_STANDARD is supported from CMake 3.1
 	set(CMAKE_CXX_STANDARD 11) # for std::unordered_set, std::unique_ptr
 	set(CMAKE_CXX_STANDARD_REQUIRED ON)
    endif() 

	add_subdirectory("${PROJECT_SOURCE_DIR}/lib")
	add_subdirectory("${PROJECT_SOURCE_DIR}/bin")
endif( ${UNIX} )
