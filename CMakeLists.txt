cmake_minimum_required(VERSION 3.10)
project(SStree)

# USO:
# 'make ss_tree_interface' para compilar solo el ejecutable de pruebas.
# 'make ss_tree_indexing' para compilar solo el ejecutable de indexación.

# Establecer el estándar C++ a utilizar
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

link_directories(/usr/lib/x86_64-linux-gnu/ /usr/lib/x86_64-linux-gnu/hdf5/serial/)

# Archivos para la rutina de indexación
set(INDEXING_SOURCE_FILES
        indexing.cpp
        params.h
        Point.h
        SStree.cpp
        SStree.h
)

# Archivos para la rutina de interfaz
set(INTERFACE_SOURCE_FILES
        Interface.cpp
        CortexAPI.cpp
        params.h
        SStree.cpp
        tinyfiledialogs.c
        CortexAPI.h
        Point.h
        SStree.h
        tinyfiledialogs.h
)

# Crear el ejecutable para la rutina de indexación
add_executable(ss_tree_indexing ${INDEXING_SOURCE_FILES})
target_include_directories(ss_tree_indexing PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Crear el ejecutable para la rutina de interfaz
add_executable(ss_tree_interface ${INTERFACE_SOURCE_FILES})
target_include_directories(ss_tree_interface PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

find_package(CURL REQUIRED)
find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(SFML 2.5 COMPONENTS graphics network REQUIRED)

target_link_libraries(ss_tree_indexing PRIVATE ${HDF5_CXX_LIBRARIES})
target_include_directories(ss_tree_indexing PRIVATE ${HDF5_CXX_INCLUDE_DIRS})

target_link_libraries(ss_tree_interface
        PRIVATE
        ${CURL_LIBRARIES}
        sfml-graphics
        sfml-network
        ${HDF5_CXX_LIBRARIES}
)

target_include_directories(ss_tree_interface
        PRIVATE
        ${CURL_INCLUDE_DIRS}
        ${HDF5_CXX_INCLUDE_DIRS}
)
