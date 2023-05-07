# find_package(Boost REQUIRED) 
# include_directories(${Boost_INCLUDE_DIRS})

add_executable(
        interp_detect
        src/interp_detect.cpp
)

target_include_directories(
        interp_detect
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/atifact>
)
target_link_libraries(interp_detect  SZ3 )
# target_link_libraries(interp_detect PUBLIC ${Boost_LIBRARIES})

add_executable(
        test_block_detection
        src/test_block_detection.cpp
)
target_include_directories(
        test_block_detection
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/atifact>
)
target_link_libraries(test_block_detection  SZ3)


add_executable(interp_detect_no_orig
        src/interp_detect_no_orig.cpp
        )
target_include_directories(
        interp_detect_no_orig
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/atifact>
)
target_link_libraries(interp_detect_no_orig  SZ3 )

add_executable(
        test_posterization
        src/test_posterization.cpp
)
target_include_directories(
        test_posterization
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/atifact>
)
target_link_libraries(test_posterization  SZ3)

install(TARGETS interp_detect EXPORT artifact
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/artifact)
