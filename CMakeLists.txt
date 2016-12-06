cmake_minimum_required(VERSION 3.1)

project( ngt )

set(ngt_VERSION_MAJOR 1 )
set(ngt_VERSION_MINOR 0 )
set(ngt_VERSION_PATCH 0 )

set( ngt_VERSION ${ngt_VERSION_MAJOR}.${ngt_VERSION_MINOR}.${ngt_VERSION_PATCH} )
set( ngt_SOVERSION ${ngt_VERSION_MAJOR} )

if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if( ${UNIX} )
	option(WALL "enable all warnings" ON)
	if( ${WALL} )
		add_compile_options(-Wall)
	endif()

	# CMAKE_CXX_STANDARD is supported from CMake 3.1
	set(CMAKE_CXX_STANDARD 11) # for std::unordered_set, std::unique_ptr
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	add_subdirectory("${PROJECT_SOURCE_DIR}/lib")
	add_subdirectory("${PROJECT_SOURCE_DIR}/bin")
endif( ${UNIX} )
