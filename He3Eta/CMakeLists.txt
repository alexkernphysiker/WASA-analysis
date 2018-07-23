add_executable(beam_momenta_offset beam_momenta_offset.cpp he3eta.cpp)
target_link_libraries(beam_momenta_offset ${ROOT_LIBRARIES}  experiment kinematics parameters)
add_executable(forward-cuts forward_cuts.cpp he3eta.cpp)
target_link_libraries(forward-cuts ${ROOT_LIBRARIES}  experiment kinematics parameters)

add_executable(luminosity_f luminosity.cpp he3eta.cpp)
target_link_libraries(luminosity_f ${ROOT_LIBRARIES} FitGen  experiment kinematics parameters)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/he3eta.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
