cmake_minimum_required(VERSION 3.25.0)
project(Project_1_EL7401 
        VERSION 0.1.0 
        DESCRIPTION "EL7041 Project 1"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # Bring in the fmt library
# include(FetchContent)
# FetchContent_Declare(
#     fmt
#     GIT_REPOSITORY https://github.com/fmtlib/fmt.github
#     GIT_TAG 10.1.1
# )

# FetchContent_MakeAvailable(fmt)

add_executable(Project_1_EL7401 main.cpp 
                src/ber.cpp
                src/data_manager.cpp
                src/eye_diagram.cpp
                src/fft.cpp
                src/impulses.cpp
                src/questions.cpp)
# target_link_libraries(Project_1_EL7401 PRIVATE fmt::fmt)

target_include_directories(Project_1_EL7401 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

