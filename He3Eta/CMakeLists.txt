add_executable(reconstruction reconstruction.cpp)
target_link_libraries(reconstruction ${ROOT_LIBRARIES} reconstructionfit experiment kinematics)
add_executable(beam_momenta_offset beam_momenta_offset.cpp he3eta.cpp)
target_link_libraries(beam_momenta_offset ${ROOT_LIBRARIES} reconstructionfit experiment kinematics)

add_executable(forward forward.cpp he3eta.cpp)
target_link_libraries(forward ${ROOT_LIBRARIES} reconstructionfit experiment kinematics)
add_executable(central central.cpp he3eta.cpp)
target_link_libraries(central ${ROOT_LIBRARIES} reconstructionfit experiment kinematics)
