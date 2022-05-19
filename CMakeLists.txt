cmake_minimum_required(VERSION 3.10)

project(mdbg)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")

add_executable(sim src/sim/main.cpp)
add_executable(mdbg 
  src/main.cpp 
  src/mdbg/minimizers.cpp
  src/mdbg/construction.cpp
  src/mdbg/opt.cpp
  src/mdbg/io.cpp)

add_executable(resolve
  src/resolve/main.cpp)

include(FetchContent)

FetchContent_Declare(
  fast
  GIT_REPOSITORY https://github.com/tbrekalo/fast.git
  GIT_TAG v0.1.0) 

FetchContent_MakeAvailable(fast)

find_package(Threads REQUIRED)

target_link_libraries(sim PRIVATE fast biosoup::biosoup)
target_link_libraries(mdbg PRIVATE fast biosoup::biosoup Threads::Threads)

target_include_directories(sim PRIVATE "include" "vendor/biosoup/include")

target_include_directories(mdbg PRIVATE 
  "include" 
  "vendor/ntHash"
  "vendor/biosoup/include" 
  "vendor/cxxopts/include"
  "vendor/robin-map/include"
  "vendor/thread-pool")

target_include_directories(resolve PRIVATE
  "include"
  "vendor/cxxopts/include"
  "vendor/robin-map/include")

## tests

add_executable(test
  test/nucleic_acid_iter.cpp
  test/custom_hash.cpp)

find_package(Catch2 2 REQUIRED)

target_link_libraries(test PRIVATE Catch2::Catch2WithMain fast biosoup::biosoup)
target_include_directories(test PRIVATE 
  "include" "vendor/ntHash" "vendor/biosoup/include")
