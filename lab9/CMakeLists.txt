set(SOURCES
    lab9.cpp
)

add_executable(lab9 ${SOURCES})
target_link_libraries(lab9 PRIVATE MPI::MPI_CXX)
