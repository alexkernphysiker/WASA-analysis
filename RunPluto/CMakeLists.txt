cmake_minimum_required(VERSION 2.6)
project(runMonteCarlo)
add_definitions(-DUSE_RANDOM_DEVICE)
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
link_directories($ENV{PLUTOSYS} $ENV{ROOTSYS}/lib)
set(CMAKE_CXX_FLAGS "-std=gnu++11")
add_executable(runhe3eta he3eta.cpp)
target_link_libraries(runhe3eta Pluto Core Gui GQt Gpad Hist RIO Proof TreePlayer EG m dl Physics)
