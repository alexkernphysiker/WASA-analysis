cmake_minimum_required(VERSION 2.6)
project(WASAatCOSY_etamesic_3He_search)
set(CMAKE_CXX_FLAGS "-std=gnu++11")
#foreign software
link_directories($ENV{PLUTOSYS} $ENV{ROOTSYS}/lib)
set(plutolibs Pluto Core Gui GQt Gpad Hist RIO Proof TreePlayer EG m dl Physics)
set(rootlibs Core Gui GQt Gpad Hist RIO Proof TreePlayer EG m dl Physics)
include_directories(
	$ENV{ROOTSYS}/include
	$ENV{PLUTOSYS}/src
	$ENV{PLUTOSYS}/plugins/brems
	$ENV{PLUTOSYS}/plugins/dalitz_mod
	$ENV{PLUTOSYS}/plugins/elementary
	$ENV{PLUTOSYS}/plugins/eta_decays
	$ENV{PLUTOSYS}/plugins/fairroot
	$ENV{PLUTOSYS}/plugins/hades
	$ENV{PLUTOSYS}/plugins/nucleus_fermi
	$ENV{PLUTOSYS}/plugins/pdg_unigen
	$ENV{PLUTOSYS}/plugins/pion_beam
	$ENV{PLUTOSYS}/plugins/scatter_mod
	$ENV{PLUTOSYS}/plugins/strangeness
	$ENV{PLUTOSYS}/plugins/w_decays
)
#libraries
add_subdirectory(lib/FitGen)
include_directories(lib/FitGen/include/)
include_directories(lib/FitGen/math_h/include)
add_subdirectory(lib/Integration)
include_directories(lib/Integration)
add_subdirectory(lib/Reconstruction)
include_directories(lib/Reconstruction)
add_subdirectory(lib/Reaction)
include_directories(lib/Reaction)
set(libs integration FitGen)

#applications
add_executable(run_plutto pluto.cpp)
target_link_libraries(run_plutto ${plutolibs} ${libs})
add_executable(reconstruction_he3 reconstruction_he3.cpp)
target_link_libraries(reconstruction_he3 ${rootlibs} ${libs} reconstruction)
add_executable(luminosity_he3eta luminosity_he3eta.cpp)
target_link_libraries(luminosity_he3eta ${rootlibs} ${libs} reaction)
