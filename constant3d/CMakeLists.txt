cmake_minimum_required (VERSION 3.20)

set (EXENAME "3DVer_constant")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#### extern
list(PREPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


add_executable (${EXENAME} 
	"main.cpp"
    "offset_crust.cpp"  "utils.cpp")

#### link libraries
find_package(CGAL REQUIRED)
target_link_libraries(${EXENAME} PRIVATE CGAL::CGAL)


find_package(TBB CONFIG REQUIRED)
if(TBB_FOUND AND NOT TARGET CGAL::TBB_support)
    if(NOT TARGET Threads::Threads)
        find_package(Threads REQUIRED)
    endif()
    add_library(CGAL::TBB_support INTERFACE IMPORTED)
    set_target_properties(CGAL::TBB_support PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_TBB;NOMINMAX"
        INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "TBB::tbb;TBB::tbbmalloc;Threads::Threads"
    )
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25)
        set_target_properties(CGAL::TBB_support PROPERTIES SYSTEM TRUE)
    endif()

    target_link_libraries(${EXENAME} PRIVATE CGAL::TBB_support)
endif()

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${EXENAME} PRIVATE Eigen3::Eigen)


find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

include(libigl)
igl_include(glfw)
target_link_libraries(${EXENAME} PRIVATE igl::glfw)

include(cli11)
target_link_libraries(${EXENAME} PRIVATE CLI11::CLI11)

#### extern libs
target_link_libraries(${EXENAME} PRIVATE PQP)
