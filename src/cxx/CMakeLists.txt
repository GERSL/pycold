option(BUILD_SCCD "Enable sccd c library" TRUE)
if (BUILD_SCCD)

  set(THREADS_PREFER_PTHREAD_FLAG ON)

  # Find required external libraries
  find_package(Threads REQUIRED)
  find_package(ZLIB REQUIRED)
  find_package(GSL REQUIRED)

  # enable_language(Fortran)
  # Corresponds to FFLAGS
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fPIC")

  # Hacks:
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  
  add_library(GLMnet STATIC "GLMnet.f")
  # set_target_properties(GLMnet PROPERTIES LINKER_LANGUAGE Fortran)

  set(SCCD_MODULE_NAME "sccd")
  # Add other C sources
  list(APPEND sccd_sources "cold.c" "input.c" "2d_array.c" "utilities.c" "misc.c" "multirobust.c" "output.c" "s_ccd.c" "KFAS.c" "lbfgs.c" "distribution_math.c")
  ## Create C++ library. Specify include dirs and link libs as normal
  #add_library(${SCCD_MODULE_NAME} SHARED ${sccd_sources})

  # HACK: statically link the library to the cython module because
  # I'm having trouble making the shared library work.
  add_library(${SCCD_MODULE_NAME} STATIC ${sccd_sources})

  #target_include_directories(
  #  ${SCCD_MODULE_NAME}
  #  PUBLIC
  #      ${NumPy_INCLUDE_DIRS}
  #      ${PYTHON_INCLUDE_DIRS}
  #      ${CMAKE_CURRENT_SOURCE_DIR}
  #)

  #message(STATUS "ZLIB::ZLIB = ${ZLIB::ZLIB}")
  #message(STATUS "Threads::Threads = ${Threads::Threads}")
  message(STATUS "Threads = ${Threads}")

 # LIB = -L$(GSL_SCI_LIB) -lz -lpthread -lgsl -lgslcblas -lgfortran -lm
  target_link_libraries(
    ${SCCD_MODULE_NAME} PUBLIC
    ZLIB::ZLIB
    Threads::Threads
    GSL::gsl
    GSL::gslcblas
    GLMnet
  )

  #target_compile_definitions(${SCCD_MODULE_NAME} PUBLIC
  #  "NPY_NO_DEPRECATED_API"
  #  #"NPY_1_7_API_VERSION=0x00000007"
  #  )

  # Transform the C++ library into an importable python module
  #python_extension_module(${SCCD_MODULE_NAME})

  # Install the C++ module to the correct relative location
  # (this will be an inplace build if you use `pip install -e`)
  #file(RELATIVE_PATH _install_dest "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  # Is this right?
  #install(TARGETS ${SCCD_MODULE_NAME} LIBRARY DESTINATION "${_install_dest}../tool/python/pycold")

else()  # build standalone program that accepts csv as input, mainly for testing purpose
  set(THREADS_PREFER_PTHREAD_FLAG ON)

  # Find required external libraries
  find_package(Threads REQUIRED)
  find_package(ZLIB REQUIRED)
  find_package(GSL REQUIRED)

  enable_language(Fortran)

  set(CMAKE_BUILD_TYPE Debug)

  # Corresponds to FFLAGS
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fPIC")

  # Hacks:
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  
  add_library(GLMnet STATIC "GLMnet.f")
  set_target_properties(GLMnet PROPERTIES LINKER_LANGUAGE Fortran)

  
  set(SCCD_MODULE_NAME "cold")
  # Add other C sources
  list(APPEND sccd_sources "cold.c" "input.c" "2d_array.c" "utilities.c" "misc.c" "multirobust.c" "output.c" "s_ccd.c" "KFAS.c" "lbfgs.c" "distribution_math.c" "sccd-desktop.c")
  ## Create C++ library. Specify include dirs and link libs as normal

  add_executable(${SCCD_MODULE_NAME} ${sccd_sources})

  #target_include_directories(
  #  ${SCCD_MODULE_NAME}
  #  PUBLIC
  #      ${NumPy_INCLUDE_DIRS}
  #      ${PYTHON_INCLUDE_DIRS}
  #      ${CMAKE_CURRENT_SOURCE_DIR}
  #)

  #message(STATUS "ZLIB::ZLIB = ${ZLIB::ZLIB}")
  #message(STATUS "Threads::Threads = ${Threads::Threads}")
  message(STATUS "Threads = ${Threads}")

  # LIB = -L$(GSL_SCI_LIB) -lz -lpthread -lgsl -lgslcblas -lgfortran -lm
  target_link_libraries(
    ${SCCD_MODULE_NAME} PUBLIC
    ZLIB::ZLIB
    Threads::Threads
    GSL::gsl
    GSL::gslcblas
    GLMnet
  )
endif()

