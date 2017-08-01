add_executable(ppn-spectra ppn_spectra.cpp)
target_link_libraries(ppn-spectra ${ROOT_LIBRARIES} gpwrap FitGen experiment kinematics)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/crosssections DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
