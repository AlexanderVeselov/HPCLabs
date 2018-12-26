set(SOURCES
    backend.cpp
    frontend.py
)

add_library(course SHARED ${SOURCES})
target_link_libraries(course PRIVATE OpenMP::OpenMP_CXX)
set(CMAKE_CXX_FLAGS "-msse2 -ftree-vectorize -fopt-info-vec-missed")
