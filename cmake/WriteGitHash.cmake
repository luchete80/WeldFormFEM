execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

file(WRITE ${OUT_HEADER} "// Auto-generated file - do not edit\n")
file(APPEND ${OUT_HEADER} "#pragma once\n")
file(APPEND ${OUT_HEADER} "#define GIT_COMMIT_HASH \"${GIT_HASH}\"\n")
file(APPEND ${OUT_HEADER} "#define PROJECT_VERSION \"${VERSION}\"\n")
string(TIMESTAMP BUILD_DATE "%Y-%m-%d")
file(APPEND ${OUT_HEADER} "#define BUILD_DATE \"${BUILD_DATE}\"\n")
