cmake_minimum_required(VERSION 3.25)

project(ClapPluginCppTemplate LANGUAGES CXX VERSION 0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flax-vector-conversions")

include(cmake/Dependencies.cmake)

# Add the object file
set(GUI_OBJ src/gui_mac.o)

add_library(${PROJECT_NAME} MODULE
    src/Utils.h
    src/Plugin.h
    src/Plugin.cpp
    src/Factory.cpp
    ${GUI_OBJ} # Include the precompiled GUI object file
)

target_link_libraries(${PROJECT_NAME} PUBLIC clap clap-helpers)
target_link_libraries(ClapPluginCppTemplate PRIVATE "-framework Accelerate")

if (APPLE)
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -Wno-unused-parameter)
    # Link Cocoa framework
    find_library(COCOA_FRAMEWORK Cocoa REQUIRED)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${COCOA_FRAMEWORK})
endif()

# Packaging
if (APPLE)
    set_target_properties(${PROJECT_NAME} PROPERTIES
        BUNDLE True
        BUNDLE_EXTENSION clap
        MACOSX_BUNDLE_GUI_IDENTIFIER org.earthics.${PROJECT_NAME}
        MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME}
        MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
        MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/cmake/ClapPluginCppTemplate.plist.in
    )
endif()

