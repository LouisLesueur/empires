cmake_minimum_required(VERSION 3.0)
project( EmpiresSim  )
find_package( OpenCV REQUIRED  )
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories( ${OpenCV_INCLUDE_DIRS}  )
add_executable( Empires src/politics.cpp src/politics.cpp 
	                src/main.cpp  )
target_link_libraries( Empires ${OpenCV_LIBS}  )
target_link_libraries (Empires Eigen3::Eigen)
