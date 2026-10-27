cmake_minimum_required(VERSION 3.16)
project(jlmfem LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Variant of the generated bindings in gen/<JLMFEM_VARIANT>: seq or par.
set(JLMFEM_VARIANT "seq" CACHE STRING "Binding variant (seq or par)")

find_package(JlCxx REQUIRED)
find_package(MFEM REQUIRED CONFIG)

set(GEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gen/${JLMFEM_VARIANT})
file(GLOB GENERATED_SOURCES CONFIGURE_DEPENDS ${GEN_DIR}/src/*.cxx)

add_library(jlmfem SHARED ${GENERATED_SOURCES})
target_include_directories(jlmfem PRIVATE ${GEN_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(jlmfem PRIVATE JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl mfem)
set_target_properties(jlmfem PROPERTIES INSTALL_RPATH "$ORIGIN")

install(TARGETS jlmfem LIBRARY DESTINATION lib)
