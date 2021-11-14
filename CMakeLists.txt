cmake_minimum_required(VERSION 3.17)
project(LRTNG)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost 1.75.0 COMPONENTS program_options)

if(Boost_FOUND)
	include_directories(/usr/local/include)
	include_directories(/usr/local/include/ibex)
	include_directories(/usr/local/include/ibex/3rd)
	include_directories(/usr/local/include/eigen3)
	include_directories(${Boost_INCLUDE_DIRS}) 

	link_directories(/usr/local/lib)
	link_directories(/usr/local/lib/ibex/3rd)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -libex -lprim")

	# used to run with AddressSanitizer, a fast memory error detector
	# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -libex -lprim -fsanitize=address -O1 -fno-omit-frame-pointer -g")

	add_executable(LRTNG main.cpp)

	target_link_libraries(LRTNG ${Boost_LIBRARIES})
endif()
