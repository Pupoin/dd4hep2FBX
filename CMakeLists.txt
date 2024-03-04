cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project("dd4hep2FBX")
set(CMAKE_BUILD_TYPE "DEBUG" )


include_directories("/home/wln/DD4hep/include/" "include/")
# set(LIB_SRC dd4hep2FBXWriter.cc)
# add_library(libdd4hep2FBX ${LIB_SRC})

aux_source_directory(./src Lib)
Add_library(libdd4hep2FBX ${Lib})


set(exename dd4hep2FBX)
# include_directories(${PROJECT_SOURCE_DIR}/dd4hep2FBX)
add_executable(${exename} test.cc)
target_link_libraries(${exename} DD4hep::DDCore DD4hep::DDCond 
        DD4hep::DDAlign ROOT::Core ROOT::Geom ROOT::GenVector CLHEP::CLHEPS CLHEP::CLHEP
        libdd4hep2FBX)

# install(TARGETS ${exename}
#   # IMPORTANT: Add the bar executable to the "export-set"
#   EXPORT ${PROJECT_NAME}Targets
#   RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)