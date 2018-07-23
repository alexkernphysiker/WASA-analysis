add_executable(luminosity_c_v2 luminosity2.cpp)
target_link_libraries(luminosity_c_v2 ${ROOT_LIBRARIES}  FitGen experiment kinematics parameters)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/crosssections DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
