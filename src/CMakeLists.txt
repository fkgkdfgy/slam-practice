add_library(main_lib SHARED
        common.cpp
        CSM.cpp
        M2MPGM.cpp
        PGM.cpp
        PoseGraph.cpp
        slam_2d_v01.cpp
        #SensorFusion.cpp
        CeresOp.cpp
        Tree.cpp
        Test.cpp)
add_executable(slam_v02 main.cpp)
target_link_libraries(slam_v02 main_lib ${the_third_party_libs})