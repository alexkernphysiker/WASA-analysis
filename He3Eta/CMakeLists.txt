add_executable(beam_momenta_offset beam_momenta_offset.cpp he3eta.cpp)
target_link_libraries(beam_momenta_offset ${ROOT_LIBRARIES} gpwrap experiment kinematics)
add_executable(forward-cuts forward_cuts.cpp he3eta.cpp)
target_link_libraries(forward-cuts ${ROOT_LIBRARIES} gpwrap experiment kinematics)

add_executable(forward forward.cpp he3eta.cpp)
target_link_libraries(forward ${ROOT_LIBRARIES} FitGen gpwrap experiment kinematics)
add_executable(forward2 forward2.cpp he3eta.cpp)
target_link_libraries(forward2 ${ROOT_LIBRARIES} FitGen gpwrap experiment kinematics)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/he3eta.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
