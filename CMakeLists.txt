cmake_minimum_required(VERSION 3.20)

project(HackathonGIG2026
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(CGAL REQUIRED)
find_package(Open3D REQUIRED)
find_package(Eigen3 REQUIRED)

# Tự động lấy tất cả file .cpp trong src/
file(GLOB_RECURSE APP_SOURCES CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

# Có thể thêm header để IDE hiện đẹp hơn
file(GLOB_RECURSE APP_HEADERS CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
)

add_executable(my_app
    ${APP_SOURCES}
    ${APP_HEADERS}
)

target_include_directories(my_app
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(my_app
    PRIVATE
        CGAL::CGAL
        Open3D::Open3D 
        Eigen3::Eigen
)

if(MSVC)
    target_compile_options(my_app PRIVATE /W4 /permissive-)
else()
    target_compile_options(my_app PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(WIN32)
    target_compile_definitions(my_app PRIVATE PLATFORM_WINDOWS)
elseif(APPLE)
    target_compile_definitions(my_app PRIVATE PLATFORM_MACOS)
elseif(UNIX)
    target_compile_definitions(my_app PRIVATE PLATFORM_LINUX)
endif()

set_target_properties(my_app PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
)