set(SOURCES
    lab6.cpp
)

add_executable(lab6 ${SOURCES})
target_link_libraries(lab6 PRIVATE OpenMP::OpenMP_CXX)
