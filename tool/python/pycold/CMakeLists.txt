option(BUILD_PYCOLD "Enable _pycold_cython module" TRUE)
if (BUILD_PYCOLD)

  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)

  set(cython_source "_pycold_cython.pyx")
  set(PYCOLD_MODULE_NAME "_pycold_cython")

  # Translate Cython into C/C++
  add_cython_target(${PYCOLD_MODULE_NAME} "${cython_source}" C OUTPUT_VAR sources)

  # Add other C sources
  list(APPEND sources )

  # Create C++ library. Specify include dirs and link libs as normal
  add_library(${PYCOLD_MODULE_NAME} MODULE ${sources})
  target_include_directories(
    ${PYCOLD_MODULE_NAME}
    PUBLIC
        ${NumPy_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_link_libraries(${PYCOLD_MODULE_NAME} ${SCCD_MODULE_NAME})

  target_compile_definitions(${PYCOLD_MODULE_NAME} PUBLIC
    "NPY_NO_DEPRECATED_API"
    #"NPY_1_7_API_VERSION=0x00000007"
    )

  # Transform the C++ library into an importable python module
  python_extension_module(${PYCOLD_MODULE_NAME})

  # Install the C++ module to the correct relative location
  # (this will be an inplace build if you use `pip install -e`)
  file(RELATIVE_PATH _install_dest "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  install(TARGETS ${PYCOLD_MODULE_NAME} LIBRARY DESTINATION "${_install_dest}")
endif()

