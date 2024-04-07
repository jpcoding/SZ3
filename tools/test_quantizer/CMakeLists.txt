add_executable(
        test_quantizer
        test_quantizer.cpp
)
target_link_libraries(test_quantizer PRIVATE SZ3)

add_executable(
        test_quantize_with_pred
        test_quantize_with_pred.cpp
)
target_link_libraries(test_quantize_with_pred PRIVATE SZ3)

target_include_directories(
    test_quantizer
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SZ3c>
)

install(TARGETS test_quantizer  EXPORT test_quantizer
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test_quantizer)
