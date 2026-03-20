cmake_minimum_required(VERSION 3.25)

project(Spectrogram LANGUAGES C CXX OBJC VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flax-vector-conversions")

include(cmake/Dependencies.cmake)

add_library(${PROJECT_NAME} MODULE
    src/Utils.h
    src/Plugin.h
    src/Plugin.cpp
    src/Factory.cpp
    src/gui_mac.m
)

target_link_libraries(${PROJECT_NAME} PUBLIC clap clap-helpers)
target_link_libraries(Spectrogram PRIVATE "-framework Accelerate")

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
        MACOSX_BUNDLE_GUI_IDENTIFIER org.earthics.Spectrogram
        MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME}
        MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
        MACOSX_BUNDLE_COPYRIGHT "Copyright © 2025 Alexander Henderson"
        MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/cmake/Spectrogram.plist.in
    )

    # Install to the user CLAP folder so hosts (e.g. Bitwig) can find the plugin.
    # Run with: cmake --install build
    install(DIRECTORY "${CMAKE_BINARY_DIR}/${PROJECT_NAME}.clap"
        DESTINATION "$ENV{HOME}/Library/Audio/Plug-Ins/CLAP"
        USE_SOURCE_PERMISSIONS
    )
endif()

