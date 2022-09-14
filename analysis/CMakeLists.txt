project(analysis)

find_package(ROOT REQUIRED)
find_package(HepMC3)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${FASTJET_DIR}/include
  ${HEPMC_DIR}/include
)
file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cc
)
file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

include(${ROOT_USE_FILE})
add_executable(analysis analysis.cc ${sources} ${headers})
add_executable(process process.cc ${sources} ${headers})
add_executable(hep hep.cc ${sources} ${headers})
add_executable(onlyptc onlyptc.cc ${sources} ${headers})
add_executable(buf buf.cc ${sources} ${headers})
add_executable(label label.cc ${sources} ${headers})
add_executable(image image.cc ${sources} ${headers})
add_executable(integrity integrity.cc ${sources} ${headers})
add_executable(recotest recotest.cc ${sources} ${headers})
#add_executable(JER JER.cc ${sources} ${headers})
#add_executable(calib calib.cc ${sources} ${headers})
target_link_libraries(
  analysis
  ${HEPMC3_LIBRARIES}
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
)
target_link_libraries(
  process
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
target_link_libraries(
  hep
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
target_link_libraries(
  onlyptc
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
target_link_libraries(
  label
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
target_link_libraries(
  buf
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
target_link_libraries(
  image
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
target_link_libraries(
  integrity
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
target_link_libraries(
  recotest
  ${HEPMC_DIR}/lib64/libHepMC3.so
  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS}
  DetComponents
  ddDRsegmentation
  stdc++fs
)
#target_link_libraries(
#  JER
#  ${HEPMC_DIR}/lib64/libHepMC3.so
#  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
#  rootIO
#  ${ROOT_LIBRARIES}
#  ${CMAKE_DL_LIBS}
#)
#target_link_libraries(
#  calib
#  ${HEPMC_DIR}/lib64/libHepMC3.so
#  ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
#  rootIO
#  ${ROOT_LIBRARIES}
#  ${CMAKE_DL_LIBS}
#)

set(
  analysis_SCRIPTS
)

foreach(_script ${analysis_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()

install(TARGETS process hep label buf image integrity recotest DESTINATION bin) # JER calib
