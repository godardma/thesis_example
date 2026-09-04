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

  find_package(CGAL REQUIRED)

  add_executable(b26_cora_inner b26_cora_inner.cpp)
  target_compile_options(b26_cora_inner PUBLIC ${CODAC_CXX_FLAGS})
  target_include_directories(b26_cora_inner SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} ${CODAC_CAPD_INCLUDE_DIR})
  target_link_libraries(b26_cora_inner PRIVATE
    ${CODAC_LIBRARIES}
    ${CODAC_CAPD_LIBRARY}
    capd::capd
    Ibex::ibex
    CGAL::CGAL
  )