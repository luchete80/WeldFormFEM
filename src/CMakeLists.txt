
add_subdirectory(common)
include_directories(../include/common)
if (BUILD_EXPLICIT)
  add_subdirectory(explicit)
endif()

if (BUILD_IMPLICIT)
  add_subdirectory(implicit)
endif()

###################################### COMMIT
set(GIT_COMMIT_HEADER "${CMAKE_BINARY_DIR}/generated/git_commit.h")

# Create directory for generated headers
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated")

# Custom command to generate git_commit.h with the short commit hash
set(GIT_COMMIT_HEADER "${CMAKE_BINARY_DIR}/generated/git_commit.h")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated")

message(STATUS "Version: ${PROJECT_VERSION}")

#add_custom_target(update_git_commit ALL DEPENDS ${GIT_COMMIT_HEADER})
add_custom_target(update_git_commit ALL
  COMMAND ${CMAKE_COMMAND} -DOUT_HEADER=${GIT_COMMIT_HEADER} -DVERSION=${PROJECT_VERSION} -P ${CMAKE_SOURCE_DIR}/cmake/WriteGitHash.cmake
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Forcing regeneration of git_commit.h"
)
add_dependencies(WeldFormFEM update_git_commit)
target_include_directories(WeldFormFEM PRIVATE "${CMAKE_BINARY_DIR}/generated")

#############################################################
