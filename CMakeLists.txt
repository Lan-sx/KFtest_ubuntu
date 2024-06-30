CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

project(Kalv1)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR}/install)
include(CTest)

find_package(ROOT COMPONENTS RIO Tree)

##############################################################################
# Package: KalTest
#    Desc: import from ILCSoft
##############################################################################
set( DICT_CINT_DEFINITIONS "HANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
set( DICT_INPUT_DIRS geomlib kallib kaltracklib utils)
set( DICT_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})

foreach( DICT ${DICT_INPUT_DIRS} )
  list( APPEND DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT} )
  aux_source_directory( src/${DICT} LIB_SOURCES )
endforeach()

#######################message 1 
#foreach(testDir ${DICT_INPUT_DIRS})
#    message("=========> $ ${testDir} ")
#endforeach()

#message("====== ${ROOT_rootcling_CMD}")

foreach( DICT ${DICT_INPUT_DIRS} )
  file( GLOB DICT_HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT}/*.h" )
  file( GLOB LINKDEFS "${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT}/LinkDef.h" )

  if( LINKDEFS )
    list( REMOVE_ITEM DICT_HEADER_FILES "${LINKDEFS}" )
  endif()
    
  add_custom_command(OUTPUT ${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm
                     COMMAND mkdir -p ${DICT_OUTPUT_DIR}
                     COMMAND ${ROOT_rootcling_CMD}
                     ARGS -f ${DICT_OUTPUT_DIR}/${DICT}.cxx -s ${DICT_OUTPUT_DIR}/${DICT} -inlineInputHeader
                     "$<$<BOOL:$<JOIN:${DICT_CINT_DEFINITIONS},>>:-D$<JOIN:${DICT_CINT_DEFINITIONS},;-D>>"
                     "$<$<BOOL:$<JOIN:${DICT_INCLUDE_DIRS},>>:-I$<JOIN:${DICT_INCLUDE_DIRS},;-I>>"
                     "$<JOIN:${DICT_HEADER_FILES},;>" "$<JOIN:${LINKDEFS},;>"
                     DEPENDS ${DICT_HEADER_FILES} ${LINKDEFS}
                     COMMAND_EXPAND_LISTS
                     )
  add_custom_target( ${DICT} DEPENDS ${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm )

  set_source_files_properties(${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm
                              PROPERTIES
                              GENERATED TRUE
                              COMPILE_FLAGS "-Wno-unused-function -Wno-overlength-strings"
			      )

  install(FILES ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm DESTINATION lib)
  install(FILES ${DICT_HEADER_FILES} DESTINATION include)

  #message("========> ${DICT_HEADER_FILES}")
  list( APPEND LIB_SOURCES ${DICT_OUTPUT_DIR}/${DICT}.cxx )
endforeach()


#
include_directories( ${DICT_INCLUDE_DIRS} )
include_directories( ${ROOT_INCLUDE_DIRS} )
include_directories( ${CMAKE_SOURCE_DIR} )
#
set(KalTestLib_srcs ${LIB_SOURCES})
#

#message("==============> ${KalTestLib_srcs}")

add_library(KalTestLib SHARED ${KalTestLib_srcs})
target_link_libraries(KalTestLib ${ROOT_LIBRARIES})

set(INSTALL_DIR ${CMAKE_CURRENT_LIST_DIR}/install/lib)

#gaudi_add_library(KalTestLib 
#                  SOURCES ${KalTestLib_srcs}
#		  LINK Gaudi::GaudiKernel
#                       ${ROOT_LIBRARIES}
#)

#
target_include_directories(KalTestLib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

message("================= $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include")
message("================= ${INSTALL_DIR}")
#
install(TARGETS KalTestLib
        LIBRARY DESTINATION ${INSTALL_DIR}
        ARCHIVE DESTINATION ${INSTALL_DIR})

#  EXPORT CEPCSWTargets
#  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
#  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
#  COMPONENT dev)



