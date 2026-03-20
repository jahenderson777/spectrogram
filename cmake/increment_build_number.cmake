# Read current build number, increment, write back, and generate header
file(READ "${BUILD_NUMBER_FILE}" BUILD_NUMBER)
string(STRIP "${BUILD_NUMBER}" BUILD_NUMBER)
math(EXPR BUILD_NUMBER "${BUILD_NUMBER} + 1")
file(WRITE "${BUILD_NUMBER_FILE}" "${BUILD_NUMBER}\n")
file(WRITE "${OUTPUT_HEADER}" "#pragma once\n#define BUILD_NUMBER ${BUILD_NUMBER}\n")
message(STATUS "Build number: ${BUILD_NUMBER}")
