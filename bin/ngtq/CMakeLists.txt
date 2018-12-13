if( ${UNIX} )
	include_directories("${PROJECT_BINARY_DIR}/lib")
        include_directories("${PROJECT_SOURCE_DIR}/lib")
        link_directories("${PROJECT_SOURCE_DIR}/lib/NGT")

	add_executable(ngtq_exe ngtq.cpp)
	add_dependencies(ngtq_exe ngt)
        set_target_properties(ngtq_exe PROPERTIES OUTPUT_NAME ngtq) 
	target_link_libraries(ngtq_exe ngt pthread)

	install(TARGETS ngtq_exe RUNTIME DESTINATION bin)

endif()
