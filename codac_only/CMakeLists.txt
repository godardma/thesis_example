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

add_executable(19_rob_arm_sing 19_rob_arm_sing.cpp)
target_compile_options(19_rob_arm_sing PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(19_rob_arm_sing SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(19_rob_arm_sing PUBLIC ${CODAC_LIBRARIES})

add_executable(20_rob_arm 20_rob_arm.cpp)
target_compile_options(20_rob_arm PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(20_rob_arm SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(20_rob_arm PUBLIC ${CODAC_LIBRARIES})

add_executable(21_rob_arm_unvalidated 21_rob_arm_unvalidated.cpp)
target_compile_options(21_rob_arm_unvalidated PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(21_rob_arm_unvalidated SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(21_rob_arm_unvalidated PUBLIC ${CODAC_LIBRARIES})

add_executable(23_box_chains 23_box_chains.cpp)
target_compile_options(23_box_chains PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(23_box_chains SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(23_box_chains PUBLIC ${CODAC_LIBRARIES})

add_executable(24_parallelepiped_chains 24_parallelepiped_chains.cpp)
target_compile_options(24_parallelepiped_chains PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(24_parallelepiped_chains SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(24_parallelepiped_chains PUBLIC ${CODAC_LIBRARIES})

add_executable(25_validation 25_validation.cpp)
target_compile_options(25_validation PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(25_validation SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(25_validation PUBLIC ${CODAC_LIBRARIES})