cmake_minimum_required(VERSION 2.6)
option(test "Build all tests." OFF)
project(WASAatCOSY_etamesic_3He_search)
set(CMAKE_CXX_FLAGS "-std=gnu++11")
option(debug "Compile in debug mode" OFF)
if(debug)
  set(CMAKE_BUILD_TYPE Debug)
endif()
#root
if(NOT $ENV{ROOTSYS} STREQUAL "")
	include_directories($ENV{ROOTSYS}/include/)
	link_directories($ENV{ROOTSYS}/lib/)
	set(ROOT_LIBRARIES Core Gui Gpad Graf Hist RIO Physics Tree)
	#libraries
	add_subdirectory(lib/FitGen)
	include_directories(lib/FitGen/include/)
	add_subdirectory(lib/Integration)
	include_directories(lib/Integration)
	add_subdirectory(lib/Reconstruction)
	include_directories(lib/Reconstruction)
	add_subdirectory(lib/Reaction)
	include_directories(lib/Reaction)
	set(libs FitGen integration reconstruction reaction)
	#applications
	add_executable(reconstruction_he3 reconstruction_he3.cpp)
	target_link_libraries(reconstruction_he3 ${ROOT_LIBRARIES} ${libs})
	add_executable(luminosity_he3eta_fit luminosity_he3eta_fit.cpp)
	target_link_libraries(luminosity_he3eta_fit ${ROOT_LIBRARIES} ${libs})
	add_executable(luminosity_he3eta_count luminosity_he3eta_eventcount.cpp)
	target_link_libraries(luminosity_he3eta_count ${ROOT_LIBRARIES} ${libs})
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
	if(test)
		add_subdirectory(tests)
	endif()
endif()