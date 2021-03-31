project(ddDRcalo)

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/compact DESTINATION .)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/compact DESTINATION ${CMAKE_BINARY_DIR}/bin)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/compact DESTINATION ${CMAKE_BINARY_DIR}/DRsim)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/compact DESTINATION ${CMAKE_BINARY_DIR}/Reco)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/compact DESTINATION ${CMAKE_BINARY_DIR}/analysis)
dd4hep_add_plugin(ddDRcalo SOURCES src/*.cpp
  USES DD4hep::DDCore DD4hep::DDCond ROOT::Core ROOT::Geom ROOT::GenVector ROOT::MathCore ddDRsegmentation)
install(TARGETS ddDRcalo LIBRARY DESTINATION lib)
install(DIRECTORY compact DESTINATION bin )
dd4hep_configure_scripts( ddDRcalo DEFAULT_SETUP WITH_TESTS)
