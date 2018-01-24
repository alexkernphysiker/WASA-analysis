add_executable(he3eta-distrs he3eta.cpp)
add_executable(pd-convert pd.cpp)
target_link_libraries(pd-convert kinematics FitGen)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pd DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
