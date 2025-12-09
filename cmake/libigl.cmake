if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG 5067c8b7eb85343af8f53cdeb0bfd1957ad59e04
)
FetchContent_MakeAvailable(libigl)
