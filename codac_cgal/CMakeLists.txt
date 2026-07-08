# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

cmake_minimum_required(VERSION 3.0.2)
project(peibos_example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Adding Codac

# In case you installed Codac in a local directory, you need 
# to specify its path with the CMAKE_PREFIX_PATH option.
# set(CMAKE_PREFIX_PATH "~/codac/build_install")

  find_package(CODAC REQUIRED)
  message(STATUS "Found Codac version ${CODAC_VERSION}")

# Initializing IBEX

  ibex_init_common()

# Compilation

if(FAST_RELEASE)
  add_compile_definitions(FAST_RELEASE)
  message(STATUS "You are running Codac in fast release mode. (option -DCMAKE_BUILD_TYPE=Release is required)")
endif()

# Adding PEIBOS

find_package(CGAL REQUIRED)

# add_executable(08_triangulate_parallelepiped_circle 08_triangulate_parallelepiped_circle.cpp)
# target_compile_options(08_triangulate_parallelepiped_circle PUBLIC ${CODAC_CXX_FLAGS})
# target_include_directories(08_triangulate_parallelepiped_circle SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
# target_link_libraries(08_triangulate_parallelepiped_circle PUBLIC ${CODAC_LIBRARIES} CGAL::CGAL)

# add_executable(13_peibos_triangulate 13_peibos_triangulate.cpp)
# target_compile_options(13_peibos_triangulate PUBLIC ${CODAC_CXX_FLAGS})
# target_include_directories(13_peibos_triangulate SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
# target_link_libraries(13_peibos_triangulate PUBLIC ${CODAC_LIBRARIES} CGAL::CGAL)

add_executable(25_validation 25_validation.cpp)
target_compile_options(25_validation PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(25_validation SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(25_validation PUBLIC ${CODAC_LIBRARIES} CGAL::CGAL)
