cmake_minimum_required(VERSION 3.0.2)
project(astar_path_planner)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  nav_msgs
  visualization_msgs
  geometry_msgs
)

find_package(Eigen3 REQUIRED)
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)

message("?${osqp_INCLUDE_DIRS}")

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES astar_path_planner
  CATKIN_DEPENDS roscpp nav_msgs visualization_msgs geometry_msgs
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  /usr/local/include/osqp
  /usr/local/include/OsqpEigen
)

add_executable(astar_planner src/astar_planner.cpp src/Trajectory.cpp)

target_link_libraries(astar_planner
  ${catkin_LIBRARIES}
  /usr/local/lib
  OsqpEigen::OsqpEigen
)