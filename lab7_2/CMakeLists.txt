set(SOURCES
    lab7_2.cpp
)

add_executable(lab7_2 ${SOURCES})
target_link_libraries(lab7_2 PRIVATE OpenMP::OpenMP_CXX)
