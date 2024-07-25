cmake_minimum_required(VERSION 3.10)
project(PhotoFusion)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 包含头文件目录
include_directories(../include)

# 设置LibRaw路径
set(LIBRAW_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../third-party/LibRaw")
add_subdirectory(${LIBRAW_PATH}-cmake ${CMAKE_CURRENT_BINARY_DIR}/LibRaw-cmake-build)

# 添加优化编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ffast-math -g")

# 添加静态库目标
add_library(PhotoFusionAPI STATIC ../src/PhotoFusionAPI.cpp ../src/Common.cpp ../src/ImageIO.cpp ../src/RawProcessor.cpp)

# 链接JPEG库
target_link_libraries(PhotoFusionAPI ${ly_lib_shared})

# 指定可执行文件和源文件
add_executable(simple_raw2tiff ./simple_raw2tiff.cpp)

# 链接静态库到可执行文件
target_link_libraries(simple_raw2tiff PRIVATE PhotoFusionAPI libraw::libraw)
