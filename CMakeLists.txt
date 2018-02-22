cmake_minimum_required(VERSION 2.8)

project( ngt )

set(ngt_VERSION_MAJOR 1 )
set(ngt_VERSION_MINOR 3 )
set(ngt_VERSION_PATCH 2 )

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

	set(BUILD_DATE_OPTION "-DBUILD_DATE=\"\\\"`date +'%Y/%m/%d %H:%M:%S'`\\\"\"")
	set(GIT_HASH_OPTION "-DGIT_HASH=\"\\\"`git log -1 --format='%H'`\\\"\"")
	set(GIT_DATE_OPTION "-DGIT_DATE=\"\\\"`git log -1 --format='%cd'`\\\"\"")
	set(GIT_TAG_OPTION "-DGIT_TAG=\"\\\"`git describe --abbrev=0`\\\"\"")

	set(BASE_OPTIONS "-Wall -std=gnu++0x -lrt ${BUILD_DATE_OPTION} ${GIT_HASH_OPTION} ${GIT_DATE_OPTION} ${GIT_TAG_OPTION}")
	if( ${NGT_AVX_DISABLED} )
		message(STATUS "AVX will not be used to compute distances.")
	else()
		set(BASE_OPTIONS "${BASE_OPTIONS} -mavx")
	endif()
	set(CMAKE_CXX_FLAGS_DEBUG "-g ${BASE_OPTIONS}")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${BASE_OPTIONS}")
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
