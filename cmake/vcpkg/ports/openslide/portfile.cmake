vcpkg_from_github(
  OUT_SOURCE_PATH SOURCE_PATH
  REPO Nico-Curti/openslide
  HEAD_REF master
  REF ...
  SHA512 ...
)

vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
  FEATURES
    "python"  PYTHON_WRAP
    "java"    BUILD_JAVA
)

set(ENABLE_PYTHON OFF)
if("python" IN_LIST FEATURES)
  x_vcpkg_get_python_packages(PYTHON_VERSION "3" PACKAGES numpy cython OUT_PYTHON_VAR "PYTHON3")
  set(ENABLE_PYTHON ON)
  set(ENV{PYTHON} "${PYTHON3}")
endif()

set(BUILD_JAVA OFF)
if ("java" IN_LIST FEATURES)
  set(BUILD_JAVA ON)
endif()


vcpkg_cmake_configure(
  SOURCE_PATH "${SOURCE_PATH}"
  OPTIONS ${FEATURE_OPTIONS}
    -DINSTALL_BIN_DIR:STRING=bin
    -DINSTALL_LIB_DIR:STRING=lib
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
    -DPYTHON_Openslide:BOOL=${ENABLE_PYTHON}
    -DBUILD_JAVA:BOOL=${BUILD_JAVA}
    -DBUILD_TEST:BOOL=OFF
    -DBUILD_DOCS:BOOL=OFF
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

vcpkg_fixup_pkgconfig()
file(INSTALL "${SOURCE_PATH}/LICENSE" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)
