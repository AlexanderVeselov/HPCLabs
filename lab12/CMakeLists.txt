set(SOURCES
    lab12.cpp
)

add_executable(lab12 ${SOURCES})
target_link_libraries(lab12 PRIVATE MPI::MPI_CXX)
