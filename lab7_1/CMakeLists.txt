set(SOURCES
    lab7_1.cpp
)

add_executable(lab7_1 ${SOURCES})
target_link_libraries(lab7_1 PRIVATE OpenMP::OpenMP_CXX)
