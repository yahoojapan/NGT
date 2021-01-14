if( ${UNIX} )
	include_directories("${PROJECT_BINARY_DIR}/lib")
        include_directories("${PROJECT_SOURCE_DIR}/lib")
        link_directories("${PROJECT_SOURCE_DIR}/lib/NGT")
##-#        link_directories("${PROJECT_SOURCE_DIR}/lib/NGTQG")

	add_executable(ngtqg_exe ngtqg.cpp)
	add_dependencies(ngtqg_exe ngt)
        set_target_properties(ngtqg_exe PROPERTIES OUTPUT_NAME ngtqg) ## 名前をngtqgに
##-#	target_link_libraries(ngtqg_exe ngtqg ngt pthread)
	target_link_libraries(ngtqg_exe ngt pthread)

	install(TARGETS ngtqg_exe RUNTIME DESTINATION bin)

endif()
