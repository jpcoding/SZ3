find_package(MPI REQUIRED)

add_executable(parallel parallel.cpp)

target_link_libraries(parallel
        PRIVATE
        SZ3
        MPI::MPI_CXX
)

target_include_directories(parallel
        PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

install(TARGETS parallel RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})