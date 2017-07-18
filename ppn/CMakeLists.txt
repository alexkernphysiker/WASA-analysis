add_executable(ppn-spectra ppn_spectra.cpp)
target_link_libraries(ppn-spectra ${ROOT_LIBRARIES} gpwrap FitGen experiment kinematics)
