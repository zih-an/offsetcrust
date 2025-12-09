cmake_minimum_required (VERSION 3.20)

file(GLOB SRC_FILES *.c *.cpp)
set (LIBNAME "PQP")

add_library(${LIBNAME} STATIC  ${SRC_FILES})
find_package(CGAL REQUIRED)
target_link_libraries(${LIBNAME} PRIVATE CGAL::CGAL)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${LIBNAME} PRIVATE Eigen3::Eigen)


target_include_directories(${LIBNAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

