include(FetchContent)

FetchContent_Declare(
    clap
    GIT_REPOSITORY https://github.com/free-audio/clap.git
    GIT_TAG main
    FIND_PACKAGE_ARGS NAMES clap
)

FetchContent_Declare(
    clap-helpers
    GIT_REPOSITORY https://github.com/free-audio/clap-helpers.git
    GIT_TAG main
    FIND_PACKAGE_ARGS NAMES clap-helpers
)

# Use FetchContent_Populate + add_subdirectory(... EXCLUDE_FROM_ALL) instead of
# FetchContent_MakeAvailable so that install() rules inside the fetched
# dependencies are suppressed.  This ensures `cmake --install build` only
# installs our plugin bundle, not clap headers into /usr/local.

FetchContent_GetProperties(clap)
if(NOT clap_POPULATED)
    FetchContent_Populate(clap)
    add_subdirectory(${clap_SOURCE_DIR} ${clap_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

FetchContent_GetProperties(clap-helpers)
if(NOT clap-helpers_POPULATED)
    FetchContent_Populate(clap-helpers)
    add_subdirectory(${clap-helpers_SOURCE_DIR} ${clap-helpers_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
