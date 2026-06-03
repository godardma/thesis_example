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
