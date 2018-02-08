add_definitions(-std=c++17)

set(THREADS ON)
add_subdirectory(math_h)
add_subdirectory(FitGen)
include_directories(${MATH_H_INC})
include_directories(${FITGEN_INC})

add_subdirectory(WASA-libs)
include_directories(${WASA_LIBS_INC})
add_subdirectory(InputDataConvert)

add_executable(reconstruction reconstruction.cpp)
target_link_libraries(reconstruction ${ROOT_LIBRARIES} reconstructionfit experiment kinematics)

add_subdirectory(He3Eta)
add_subdirectory(Gammas)
add_subdirectory(ppn)
