cmake_minimum_required(VERSION 3.20)  


if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()


set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)


set(CMAKE_POLICY_VERSION_MINIMUM 3.5)


#set(CMAKE_TOOLCHAIN_FILE "E:/vcpkg/scripts/buildsystems/vcpkg.cmake")
set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_BUILD_TYPE RELEASE)

project ("offset")



## PQP
add_subdirectory("extern/PQP-model3D")
install(TARGETS PQP
	EXPORT ${PROJECT_NAME}Targets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include
)	
install(
	DIRECTORY ./include/ 
	DESTINATION include
    FILES_MATCHING 
	PATTERN "*.h" 
	PATTERN "*.hpp"
)



## MAIN
add_subdirectory ("constant3d")

add_subdirectory ("variable3d")

add_subdirectory ("examples")
add_subdirectory ("evaluation")

