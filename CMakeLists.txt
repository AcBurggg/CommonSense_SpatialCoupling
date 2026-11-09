# ------------------------------- #
# -- Tianji's Project Template -- #
# ------------------------------- #

cmake_minimum_required(VERSION 3.16...3.27)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(tianji NO_POLICY_SCOPE)

# -- Your project name here

project(doro LANGUAGES CXX C)

# -- Configure project compiler options

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS true)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(ENV{CC}     ${CMAKE_C_COMPILER})
set(ENV{CXX}    ${CMAKE_CXX_COMPILER})
enable_testing()

find_package(Threads REQUIRED)

# -- Manage Compile Options w/ ASAN flags

if(CMAKE_BUILD_TYPE MATCHES "Release")
    add_compile_options(-O3 -Werror -Wall -Wno-attributes -DNDEBUG)
    # Eigen 3.4.0's ARM NEON code trips -Wclass-memaccess (not hit on x86)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-error=class-memaccess>)
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
    add_compile_options(-O0 -ggdb -Werror -Wall -Wno-attributes -DDEBUG -DONIAK_DEBUG)
endif()

# -- Process Main Library & Unit Tests

AutoBuild(
    LIB_DIR "lib" STATIC
    BIN_DIR "bin"
)

# -- Manage Third Party Dependencies

Git(
    SITE    "https://github.com"
    USER    "gflags"
    REPO    "gflags"
    BRANCH  "v2.2.2"
    PIPELINE "CMAKE INSTALL" FLAGS "-DBUILD_SHARED_LIBS=ON" PUBLIC
)
Git(
    SITE    "https://github.com"
    USER    "google"
    REPO    "glog"
    BRANCH  "v0.7.0"
    PIPELINE "CMAKE INSTALL" PUBLIC
)
Git(
    SITE    "https://github.com"
    USER    "google"
    REPO    "googletest"
    BRANCH  "v1.14.0"
    PACK    "GTest"
    TARGET  "gtest"
    PIPELINE "CMAKE SUBDIR"
)
Git(
    SITE    "https://github.com"
    USER    "nlohmann"
    REPO    "json"
    PACK    "nlohmann_json"
    BRANCH  "v3.11.3"
    PIPELINE "CMAKE INSTALL" FLAGS "-DJSON_BuildTests=OFF" PUBLIC
)
Git(
    SITE    "https://gitlab.com"
    USER    "libeigen"
    REPO    "eigen"
    PACK    "Eigen3" 
    BRANCH  "3.4.0"
    TARGET  "Eigen"
    PIPELINE "CMAKE INSTALL" PUBLIC
)
Git(
    SITE    "https://github.com"
    USER    "wangyi-fudan"
    REPO    "wyhash"
    BRANCH  "wyhash_final4"
    PIPELINE "INCLUDE" PUBLIC
    INCLUDE_DIR "${CMAKE_BINARY_DIR}/3rd_party_download/wangyi-fudan"
)
Git(
    SITE   "https://github.com"
    USER    "Y-jiji"
    REPO   "rans-fix"
    BRANCH "master"
    PIPELINE "INCLUDE" PUBLIC
    INCLUDE_DIR "${CMAKE_BINARY_DIR}/3rd_party_download/Y-jiji"
)
Git(
    SITE  "https://github.com"
    USER    "Y-jiji"
    REPO    "bch_codec"
    BRANCH  "master"
    PIPELINE "CMAKE INSTALL" PUBLIC
)
Git(
    SITE    "https://github.com"
    USER    "Ten0"
    REPO    "updatable_priority_queue"
    BRANCH  "master"
    PIPELINE "INCLUDE" PUBLIC
    INCLUDE_DIR "${CMAKE_BINARY_DIR}/3rd_party_download/Ten0"
)

# -- Specially treat libONIAK
target_include_directories(
    ${CMAKE_PROJECT_NAME}
    PUBLIC libONIAK
)