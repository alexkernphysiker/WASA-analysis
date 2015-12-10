cmake_minimum_required(VERSION 2.6)
project(WASAatCOSY_etamesic_3He_search)
set(CMAKE_CXX_FLAGS "-std=gnu++11")
#root
foreach(mode QUIET REQUIRED)
  find_package(ROOT 5 ${mode} COMPONENTS Hist Physics RIO)
  if(ROOT_USE_FILE)
    include(${ROOT_USE_FILE})
  endif()
  if(ROOT_FOUND)
    break()
  endif()
endforeach()
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_DEFINITIONS})
link_directories(${ROOT_LIBRARY_DIRS})
#libraries
add_subdirectory(lib/FitGen)
include_directories(lib/FitGen/include/)
add_subdirectory(lib/Integration)
include_directories(lib/Integration)
add_subdirectory(lib/Reconstruction)
include_directories(lib/Reconstruction)
add_subdirectory(lib/Reaction)
include_directories(lib/Reaction)
set(libs integration FitGen)
#applications
add_executable(reconstruction_he3 reconstruction_he3.cpp)
target_link_libraries(reconstruction_he3 ${ROOT_LIBRARIES} ${libs} reconstruction)
add_executable(luminosity_he3eta luminosity_he3eta.cpp)
target_link_libraries(luminosity_he3eta ${ROOT_LIBRARIES} ${libs} reaction)
#pluto
if(NOT $ENV{PLUTOSYS} STREQUAL "")
link_directories($ENV{PLUTOSYS})
	include_directories(
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
	add_executable(run_plutto pluto.cpp)
	target_link_libraries(run_plutto Pluto ${ROOT_LIBRARIES} ${libs})
endif()
option(test "Build all tests." OFF) 
if(test)
	add_subdirectory(tests)
endif()