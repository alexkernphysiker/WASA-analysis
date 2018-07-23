add_executable(central-2g central-2g.cpp)
target_link_libraries(central-2g ${ROOT_LIBRARIES} experiment kinematics FitGen parameters)
add_executable(central-6g central-6g.cpp)
target_link_libraries(central-6g ${ROOT_LIBRARIES} experiment kinematics FitGen parameters)
