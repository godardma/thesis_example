# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

  cmake_minimum_required(VERSION 3.5)
  project(codac_example LANGUAGES CXX)

  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Adding Codac

  # In case you installed Codac in a local directory, you need 
  # to specify its path with the CMAKE_PREFIX_PATH option.
  # set(CMAKE_PREFIX_PATH "~/codac/build_install")

  find_package(CODAC REQUIRED)
  message(STATUS "Found Codac version ${CODAC_VERSION}")

# Initializating Ibex
  
  ibex_init_common()

# Compilation

  if(FAST_RELEASE)
    add_compile_definitions(FAST_RELEASE)
    message(STATUS "You are running Codac in fast release mode. (option -DCMAKE_BUILD_TYPE=Release is required)")
  endif()

  add_executable(14_pendulum 14_pendulum.cpp)
  target_compile_options(14_pendulum PUBLIC ${CODAC_CXX_FLAGS})
  target_include_directories(14_pendulum SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} ${CODAC_CAPD_INCLUDE_DIR})
  target_link_libraries(14_pendulum PRIVATE
    ${CODAC_LIBRARIES}
    ${CODAC_CAPD_LIBRARY}
    capd::capd
    Ibex::ibex
  )

  add_executable(15_lorenz 15_lorenz.cpp)
  target_compile_options(15_lorenz PUBLIC ${CODAC_CXX_FLAGS})
  target_include_directories(15_lorenz SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} ${CODAC_CAPD_INCLUDE_DIR})
  target_link_libraries(15_lorenz PRIVATE
    ${CODAC_LIBRARIES}
    ${CODAC_CAPD_LIBRARY}
    capd::capd
    Ibex::ibex
  )
  
  add_executable(16_spiral 16_spiral.cpp)
  target_compile_options(16_spiral PUBLIC ${CODAC_CXX_FLAGS})
  target_include_directories(16_spiral SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} ${CODAC_CAPD_INCLUDE_DIR})
  target_link_libraries(16_spiral PRIVATE
    ${CODAC_LIBRARIES}
    ${CODAC_CAPD_LIBRARY}
    capd::capd
    Ibex::ibex
  )
  
  add_executable(17_dubins 17_dubins.cpp)
  target_compile_options(17_dubins PUBLIC ${CODAC_CXX_FLAGS})
  target_include_directories(17_dubins SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} ${CODAC_CAPD_INCLUDE_DIR})
  target_link_libraries(17_dubins PRIVATE
    ${CODAC_LIBRARIES}
    ${CODAC_CAPD_LIBRARY}
    capd::capd
    Ibex::ibex
  )

  add_executable(18_cshape 18_cshape.cpp)
  target_compile_options(18_cshape PUBLIC ${CODAC_CXX_FLAGS})
  target_include_directories(18_cshape SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} ${CODAC_CAPD_INCLUDE_DIR})
  target_link_libraries(18_cshape PRIVATE
    ${CODAC_LIBRARIES}
    ${CODAC_CAPD_LIBRARY}
    capd::capd
    Ibex::ibex
  )

  add_executable(22_area_reach 22_area_reach.cpp)
  target_compile_options(22_area_reach PUBLIC ${CODAC_CXX_FLAGS})
  target_include_directories(22_area_reach SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} ${CODAC_CAPD_INCLUDE_DIR})
  target_link_libraries(22_area_reach PRIVATE
    ${CODAC_LIBRARIES}
    ${CODAC_CAPD_LIBRARY}
    capd::capd
    Ibex::ibex
  )
