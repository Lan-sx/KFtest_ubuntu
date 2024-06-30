cmake_minimum_required(VERSION 3.0.0)
project(KalSimple)

include(CTest)
enable_testing()

#---Generate JSON file for VIM 
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#---Gdb Debug
#SET(CMAKE_BUILD_TYPE "Debug")
#SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g3 -ggdb -gdwarf-2")
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie")

add_compile_options(-fPIC)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR}/lib)
include(CTest)

find_package(ROOT REQUIRED COMPONENTS RIO Hist Tree Net Gui Spectrum)


set(KalLibdir /mnt/g/SoftWare/CERN/WSL/ILCsoft/KalTest/install/lib)
set(KalIncludedir /mnt/g/SoftWare/CERN/WSL/ILCsoft/KalTest/install/include)

include(${ROOT_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${KalIncludedir})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


add_executable(${PROJECT_NAME} EXKalTest.cxx ${sources} ${headers})

#target_include_directories(${PROJECT_NAME} PUBLIC ${KalIncludedir} ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${KalIncludedir})

#target_link_directories(${PROJECT_NAME} PRIVATE ${KalLibdir} ${CMAKE_CURRENT_LIST_DIR})
target_link_directories(${PROJECT_NAME} PRIVATE ${KalLibdir} )

target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} KalTestLib)
#target_link_libraries(${PROJECT_NAME}  KalTestLib)


#foreach(_script ${EXAMPLEB1_SCRIPTS})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_script}
#    ${PROJECT_BINARY_DIR}/${_script}
#    COPYONLY
#    )
#endforeach()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
