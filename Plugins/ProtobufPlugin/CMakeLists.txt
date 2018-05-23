file (GLOB_RECURSE src_files "src/*.cpp" "gen/*.cc" "include/*.hpp")

## PROTOBUF DEPENDENCY SETUP BEGIN

include(ExternalProject)

ExternalProject_Add(PROTOBUF_PROJECT

  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/protobuf

  GIT_REPOSITORY "https://github.com/google/protobuf.git"
  GIT_TAG "v3.5.1"
  GIT_SUBMODULES
  
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/protobuf"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} -DCMAKE_POSITION_INDEPENDENT_CODE=ON -Dprotobuf_BUILD_TESTS=OFF <SOURCE_DIR>/cmake
  INSTALL_COMMAND ""

  TEST_COMMAND ""
  STEP_TARGETS build
)

ExternalProject_Get_Property(PROTOBUF_PROJECT BINARY_DIR)
set(_PROTOBUF_BINARY_DIR ${BINARY_DIR})

set(_PROTOBUF_LIB_LOCATION "${_PROTOBUF_BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}protobuf${CMAKE_STATIC_LIBRARY_SUFFIX}")

# this is not great, but otherwise ninja won't play nice with the build
if("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  message(STATUS "Fixing dependencies for Ninja")
  add_custom_command(OUTPUT ${_PROTOBUF_LIB_LOCATION} DEPENDS PROTOBUF_PROJECT-build)
  add_custom_target(PROTOBUF_LIB DEPENDS ${_PROTOBUF_LIB_LOCATION})
endif()

add_library(protobuf STATIC IMPORTED)
set_target_properties(protobuf PROPERTIES IMPORTED_LOCATION ${_PROTOBUF_LIB_LOCATION})
add_dependencies(protobuf PROTOBUF_PROJECT)

# set standard variables used to link and include
set(PROTOBUF_LIBRARY protobuf)
ExternalProject_Get_Property(PROTOBUF_PROJECT SOURCE_DIR)
set(PROTOBUF_INCLUDE_DIRS "${SOURCE_DIR}/src")

## PROTOBUF DEPENDENCY SETUP END



add_library (ACTSProtobufPlugin SHARED ${src_files})
target_link_libraries(ACTSProtobufPlugin PRIVATE ${PROTOBUF_LIBRARY})
target_include_directories (ACTSProtobufPlugin PRIVATE ${PROTOBUF_INCLUDE_DIRS})
target_include_directories (ACTSProtobufPlugin PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ACTSProtobufTest src/test.cpp)
target_link_libraries(ACTSProtobufPlugin PRIVATE ${PROTOBUF_LIBRARY})
target_link_libraries(ACTSProtobufTest PUBLIC ACTSProtobufPlugin)
target_include_directories (ACTSProtobufTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (ACTSProtobufTest PRIVATE ${PROTOBUF_INCLUDE_DIRS})

install (TARGETS ACTSProtobufPlugin
         EXPORT ACTSProtobufPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT ProtobufPlugin TARGETS ACTSProtobufPlugin)
