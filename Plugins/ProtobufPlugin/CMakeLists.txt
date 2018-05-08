file (GLOB_RECURSE src_files "src/*.cpp" "gen/*.cc" "include/*.hpp")


#include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(ExternalProject)


ExternalProject_Add(PROTOBUF_PROJECT

  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/protobuf

  GIT_REPOSITORY "https://github.com/google/protobuf.git"
  GIT_TAG "v3.5.1"
  GIT_SUBMODULES
  
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/protobuf"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} -Dprotobuf_BUILD_TESTS=OFF <SOURCE_DIR>/cmake
  INSTALL_COMMAND ""

  TEST_COMMAND ""
)

ExternalProject_Get_Property(PROTOBUF_PROJECT BINARY_DIR)
set(PROTOBUF_BINARY_DIR ${BINARY_DIR})
ExternalProject_Get_Property(PROTOBUF_PROJECT SOURCE_DIR)
set(PROTOBUF_INCLUDE_DIRS "${SOURCE_DIR}/src")

add_library (ACTSProtobufPlugin SHARED ${src_files})
add_dependencies(ACTSProtobufPlugin PROTOBUF_PROJECT)

link_directories(${PROTOBUF_BINARY_DIR})

target_link_libraries (ACTSProtobufPlugin PRIVATE protobuf)
target_include_directories (ACTSProtobufPlugin PRIVATE ${PROTOBUF_INCLUDE_DIRS})
target_include_directories (ACTSProtobufPlugin PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ACTSProtobufTest src/test.cpp)
target_include_directories (ACTSProtobufTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ACTSProtobufTest PRIVATE protobuf)
target_link_libraries(ACTSProtobufTest PUBLIC ACTSProtobufPlugin)


install (TARGETS ACTSProtobufPlugin
         EXPORT ACTSProtobufPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#acts_add_targets_to_cdash_project(PROJECT ProtobufPlugin TARGETS ACTSProtobufPlugin)
