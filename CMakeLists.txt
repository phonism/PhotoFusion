cmake_minimum_required(VERSION 3.10)
project(PhotoFusion)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 包含头文件目录
include_directories(include)

# 添加 LibRaw-cmake 子目录
set(LIBRAW_PATH "${CMAKE_CURRENT_SOURCE_DIR}/third-party/LibRaw")
add_subdirectory(third-party/LibRaw-cmake)

# 查找 OpenMP 包
find_package(OpenMP REQUIRED)

# 添加优化编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -funroll-loops")

# 指定可执行文件和源文件
add_executable(simple_raw2tiff samples/simple_raw2tiff.cpp src/Common.cpp src/ImageIO.cpp src/RawProcessor.cpp src/PhotoFusionAPI.cpp)

# 链接 LibRaw 和 OpenMP 库
if(OpenMP_CXX_FOUND)
    target_link_libraries(simple_raw2tiff PRIVATE libraw::libraw ${OpenMP_CXX_FLAGS})
else()
    message(FATAL_ERROR "OpenMP not found")
endif()
