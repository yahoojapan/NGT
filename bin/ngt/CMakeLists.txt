if( ${UNIX} )
	include_directories("${PROJECT_BINARY_DIR}/lib")
        include_directories("${PROJECT_SOURCE_DIR}/lib")
        link_directories("${PROJECT_SOURCE_DIR}/lib/NGT")

	add_executable(ngt_exe ngt.cpp)
	add_dependencies(ngt_exe ngt)
        set_target_properties(ngt_exe PROPERTIES OUTPUT_NAME ngt) 
 	target_link_libraries(ngt_exe ngt pthread)
	install(TARGETS ngt_exe RUNTIME DESTINATION bin)

endif()
