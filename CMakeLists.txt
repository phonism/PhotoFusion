cmake_minimum_required(VERSION 3.10)
project(PhotoFusion)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 包含头文件目录
include_directories(include)

# 添加优化编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ffast-math")

# 添加静态库目标
add_library(PhotoFusionAPI STATIC src/PhotoFusionAPI.cpp src/Common.cpp src/ImageIO.cpp src/RawProcessor.cpp)

target_link_libraries(PhotoFusionAPI PUBLIC)
