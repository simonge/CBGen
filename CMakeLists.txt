#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CBGen)
set(useROOT true)        # use true or false
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#


set(DIRACXX_HOME /home/simong/Generators/GlueXCB/Diracxx-master)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(/home/simong/geant4/FindROOT.cmake)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${DIRACXX_HOME})
set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR} ${DIRACXX_HOME}
)

link_directories( ${LINK_DIRECTORIES})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB main     ${PROJECT_SOURCE_DIR}/CBGen.cc)
file(GLOB sources  ${PROJECT_SOURCE_DIR}/src/*.cc) 
file(GLOB headers  ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(CBGen ${main} ${sources} ${headers})
add_executable(Adapt adapt.cc ${sources} ${headers})
set(LINK_LIBRARIES
 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lDirac
)


target_link_libraries(CBGen ${LINK_LIBRARIES})
target_link_libraries(Adapt ${LINK_LIBRARIES})

#----------------------------------------------------------------------------

#set(EXECUTABLE_OUTPUT_PATH ) 

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CBGen DESTINATION bin)
install(TARGETS Adapt DESTINATION bin)


