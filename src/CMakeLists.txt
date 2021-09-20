option(BUILD_SCCD "Enable sccd c library" TRUE)
if (BUILD_SCCD)

  set(THREADS_PREFER_PTHREAD_FLAG ON)

  # Find required external libraries
  find_package(Threads REQUIRED)
  find_package(ZLIB REQUIRED)
  find_package(GSL REQUIRED)
  #find_package(m REQUIRED)

  set(SCCD_MODULE_NAME "sccd")

  # Add other C sources
  list(APPEND sources "cold.c" "input.c" "2d_array.c" "utilities.c" "misc.c" "multirobust.c" "output.c" "s_ccd.c" "KFAS.c" "lbfgs.c" "distribution_math.c")

  ## Create C++ library. Specify include dirs and link libs as normal
  add_library(${SCCD_MODULE_NAME} MODULE ${sources})

  #target_include_directories(
  #  ${SCCD_MODULE_NAME}
  #  PUBLIC
  #      ${NumPy_INCLUDE_DIRS}
  #      ${PYTHON_INCLUDE_DIRS}
  #      ${CMAKE_CURRENT_SOURCE_DIR}
  #)

 # LIB = -L$(GSL_SCI_LIB) -lz -lpthread -lgsl -lgslcblas -lgfortran -lm -fopenmp
  target_link_libraries(
    ${SCCD_MODULE_NAME} PRIVATE
    ZLIB::ZLIB
    Threads::Threads
    GSL::gsl
    GSL::gslcblas
    #${M_LIBRARIES}
  )

  #target_compile_definitions(${SCCD_MODULE_NAME} PUBLIC
  #  "NPY_NO_DEPRECATED_API"
  #  #"NPY_1_7_API_VERSION=0x00000007"
  #  )

  # Transform the C++ library into an importable python module
  #python_extension_module(${SCCD_MODULE_NAME})

  # Install the C++ module to the correct relative location
  # (this will be an inplace build if you use `pip install -e`)
  file(RELATIVE_PATH _install_dest "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  # Is this right?
  install(TARGETS ${SCCD_MODULE_NAME} LIBRARY DESTINATION "${_install_dest}../tool/python/pycold")

endif()

