cmake_minimum_required(VERSION 3.10)
project(SStree)

# USO:
# 'make ss_tree_interface' para compilar solo el ejecutable de pruebas.
# 'make ss_tree_indexing' para compilar solo el ejecutable de indexación.

# Establecer el estándar C++ a utilizar
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Archivos para la rutina de indexación
set(INDEXING_SOURCE_FILES
        indexing.cpp
        inc/params.h
        inc/Point.h
        src/SStree.cpp
        inc/SStree.h
)

# Archivos para la rutina de interfaz
set(INTERFACE_SOURCE_FILES
        interface.cpp
        src/CortexAPI.cpp
        inc/params.h
        src/SStree.cpp
        src/tinyfiledialogs.c
        inc/CortexAPI.h
        inc/Point.h
        inc/SStree.h
        inc/tinyfiledialogs.h
)

# Crear el ejecutable para la rutina de indexación
add_executable(ss_tree_indexing ${INDEXING_SOURCE_FILES})
target_include_directories(ss_tree_indexing PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Crear el ejecutable para la rutina de interfaz
add_executable(ss_tree_interface ${INTERFACE_SOURCE_FILES})
target_include_directories(ss_tree_interface PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

find_package(CURL REQUIRED)
find_package(SFML 2.5 COMPONENTS graphics network REQUIRED)

target_link_libraries(ss_tree_interface
        PRIVATE
        ${CURL_LIBRARIES}
        sfml-graphics
        sfml-network
)

target_include_directories(ss_tree_interface
        PRIVATE
        ${CURL_INCLUDE_DIRS}
)
