cmake_minimum_required(VERSION 3.20)

set(LIBNAME LibCore)
set(EXENAME 3DVer_variable)

list(PREPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SRC_FILES
    utils.cpp
    offset_crust_neq.cpp
)

add_library(${LIBNAME} STATIC ${SRC_FILES})

find_package(CGAL REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)
find_package(OpenMP REQUIRED)
find_package(nanoflann CONFIG REQUIRED)
include(libigl)
include(cli11)

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
endif()


if(OpenMP_CXX_FOUND)
    target_compile_options(${LIBNAME} PUBLIC ${OpenMP_CXX_FLAGS})
    target_link_libraries(${LIBNAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${LIBNAME} PUBLIC
    CGAL::CGAL
    CGAL::TBB_support
    Eigen3::Eigen
    igl::glfw
    CLI11::CLI11
    PQP
    nanoflann::nanoflann
)

target_include_directories(${LIBNAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

add_executable(${EXENAME} main.cpp)
target_link_libraries(${EXENAME} PRIVATE ${LIBNAME})
