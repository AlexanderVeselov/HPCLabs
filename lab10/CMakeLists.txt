set(SOURCES
    lab10.cpp
)

add_executable(lab10 ${SOURCES})
target_link_libraries(lab10 PRIVATE MPI::MPI_CXX)
