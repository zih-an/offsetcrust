cmake_minimum_required (VERSION 3.20)

set (EXENAME "cruxt_examples")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable (${EXENAME} 
    "DupinCyclide.cpp"
    "Knot.cpp" 
    "MAT_recon.cpp"
    "main.cpp"
)

target_link_libraries(${EXENAME} PRIVATE LibCore)

