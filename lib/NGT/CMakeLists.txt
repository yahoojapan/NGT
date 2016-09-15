if( ${UNIX} )
	option(NGT_SHARED_MEMORY_ALLOCATOR "enable shared memory" OFF)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defines.h.in ${CMAKE_CURRENT_BINARY_DIR}/defines.h)
	include_directories("${CMAKE_CURRENT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/lib" "${PROJECT_BINARY_DIR}/lib/")

	file(GLOB NGT_SOURCES *.cpp)

	add_library(ngtstatic STATIC ${NGT_SOURCES})
	set_target_properties(ngtstatic PROPERTIES OUTPUT_NAME ngt) 
	target_link_libraries(ngtstatic)


	add_library(ngt SHARED ${NGT_SOURCES})
	set_target_properties(ngt PROPERTIES VERSION ${ngt_VERSION}) 
	set_target_properties(ngt PROPERTIES SOVERSION ${ngt_SOVERSION}) 
	add_dependencies(ngt ngtstatic) 
	target_link_libraries(ngt)

	install(TARGETS ngt LIBRARY DESTINATION lib64)
	install(TARGETS ngtstatic ARCHIVE DESTINATION lib64)
endif()
