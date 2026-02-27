set(USER-BACKMAP_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB USER-BACKMAP_SOURCES ${USER-BACKMAP_SOURCES_DIR}/*.cpp)
file(GLOB USER-BACKMAP_HEADERS ${USER-BACKMAP_SOURCES_DIR}/*.h)

foreach(MY_CXX_SRC ${USER-BACKMAP_SOURCES})
  get_filename_component(MY_SRC_NAME ${MY_CXX_SRC} NAME)
  target_sources(lammps PRIVATE ${MY_CXX_SRC})
endforeach()

foreach(MY_HEADER ${USER-BACKMAP_HEADERS})
  get_filename_component(MY_HDR_NAME ${MY_HEADER} NAME)
  target_sources(lammps PRIVATE ${MY_HEADER})
endforeach()
