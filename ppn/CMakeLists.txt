add_executable(luminosity_c luminosity.cpp)
target_link_libraries(luminosity_c ${ROOT_LIBRARIES}  FitGen experiment kinematics)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/crosssections DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
