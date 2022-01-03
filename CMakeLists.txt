cmake_minimum_required(VERSION 3.10)

project(mdbg)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined")

add_executable(sim src/sim/main.cpp)
add_executable(mdbg 
  src/main.cpp 
  src/mdbg/minimizers.cpp)

include(FetchContent)

FetchContent_Declare(
  fast
  GIT_REPOSITORY https://github.com/tbrekalo/fast.git
  GIT_TAG v0.1.0) 

FetchContent_MakeAvailable(fast)

target_link_libraries(sim PRIVATE fast)
target_link_libraries(mdbg PRIVATE fast)
target_include_directories(sim PRIVATE "include")
target_include_directories(mdbg PRIVATE "include" "vendor/ntHash")
