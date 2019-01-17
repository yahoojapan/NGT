if( ${UNIX} )
	option(NGT_SHARED_MEMORY_ALLOCATOR "enable shared memory" OFF)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defines.h.in ${CMAKE_CURRENT_BINARY_DIR}/defines.h)
	include_directories("${CMAKE_CURRENT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/lib" "${PROJECT_BINARY_DIR}/lib/")

	file(GLOB NGT_SOURCES *.cpp)
	file(GLOB HEADER_FILES *.h *.hpp)
	file(GLOB NGTQ_HEADER_FILES NGTQ/*.h NGTQ/*.hpp)

	add_library(ngtstatic STATIC ${NGT_SOURCES})
	set_target_properties(ngtstatic PROPERTIES OUTPUT_NAME ngt) 
	set_target_properties(ngtstatic PROPERTIES COMPILE_FLAGS "-fPIC")
	target_link_libraries(ngtstatic)

	add_library(ngt SHARED ${NGT_SOURCES})
	set_target_properties(ngt PROPERTIES VERSION ${ngt_VERSION}) 
	set_target_properties(ngt PROPERTIES SOVERSION ${ngt_SOVERSION}) 
	add_dependencies(ngt ngtstatic) 
	if(${APPLE})
		target_link_libraries(ngt gomp)
	else(${APPLE})
		target_link_libraries(ngt gomp rt)
	endif(${APPLE})

	install(TARGETS
		ngt
		ngtstatic
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib)
	install(FILES ${HEADER_FILES}  ${CMAKE_CURRENT_BINARY_DIR}/defines.h DESTINATION include/NGT)
	install(FILES ${NGTQ_HEADER_FILES} DESTINATION include/NGT/NGTQ)
endif()
