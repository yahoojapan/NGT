if( ${UNIX} )
	include_directories("${PROJECT_BINARY_DIR}/lib")
        include_directories("${PROJECT_SOURCE_DIR}/lib")
        link_directories("${PROJECT_SOURCE_DIR}/lib/NGT")

	add_executable(search search.cpp)
	add_dependencies(search ngt)
	target_link_libraries(search ngt pthread)
endif()
