project(analysis)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

add_executable(analysis analysis.cpp ${sources} ${headers})
add_executable(process process.cpp ${sources} ${headers})

target_link_libraries(
  analysis
  ${ROOT_LIBRARIES}
  podio::podioRootIO
  edm4dr
  edm4dr::edm4drDict
)
target_link_libraries(
  process
  ${ROOT_LIBRARIES}
  podio::podioRootIO
  edm4dr
  edm4dr::edm4drDict
)
  #DRcomponents #FIXME
  #DRsegmentation

target_include_directories(analysis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_include_directories(process PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(analysis PROPERTIES PUBLIC_HEADER "${headers}")
set_target_properties(process PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS analysis process EXPORT analysisTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)
