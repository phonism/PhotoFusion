cmake_minimum_required(VERSION 2.5)
project(PhotoFusion)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加选项以启用或禁用 OpenMP
option(USE_OPENMP "Enable OpenMP support" ON)
option(USE_FFTW "Enable FFTW support" ON)
option(USE_GCD "Enable GCD support" ON)

# 包含头文件目录
include_directories(../include)

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        add_definitions(-DUSE_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

if(USE_FFTW)
    add_definitions(-DUSE_FFTW)
    # 查找FFTW
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/.")
    find_package(FFTW REQUIRED)
    include_directories(${FFTW_INCLUDE_DIRS})
    link_directories(${FFTW_LIBRARY_DIRS})
endif()

if(USE_GCD)
    add_definitions(-DUSE_GCD)
endif()

# 添加优化编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ffast-math")

# 添加静态库目标
file(GLOB SRC_FILES ../src/*.cpp)
add_library(PhotoFusionAPI STATIC ${SRC_FILES})

# 链接JPEG库
#target_link_libraries(PhotoFusionAPI ${ly_lib_shared})

# 指定可执行文件和源文件
add_executable(simple_raw2tiff ./simple_raw2tiff.cpp)

# 指定可执行文件和源文件
add_executable(simple_stack ./simple_stack.cpp)

if(OpenMP_CXX_FOUND AND USE_FFTW)
    target_link_libraries(PhotoFusionAPI ${ly_lib_shared} ${OpenMP_CXX_FLAGS} ${FFTW_LIBRARIES} fftw3_threads fftw3)
    target_link_libraries(simple_raw2tiff PRIVATE PhotoFusionAPI ${OpenMP_CXX_FLAGS} ${FFTW_LIBRARIES} fftw3_threads fftw3)
    target_link_libraries(simple_stack PRIVATE PhotoFusionAPI ${OpenMP_CXX_FLAGS} ${FFTW_LIBRARIES} fftw3_threads fftw3)
else()
    target_link_libraries(PhotoFusionAPI ${ly_lib_shared})
    target_link_libraries(simple_raw2tiff PRIVATE PhotoFusionAPI)
    target_link_libraries(simple_stack PRIVATE PhotoFusionAPI)
endif()
