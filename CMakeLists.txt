PROJECT(mysolver)
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

ADD_EXECUTABLE(main3d 
	main.cpp
        FlowGravityForcing.cpp         
        LevelSetInitialCondition.cpp   
        SetLSProperties.cpp
        GravityForcing.cpp             
        SetFluidGasSolidDensity.cpp    
        TagLSRefinementCells.cpp
        LSLocateGasInterface.cpp       
        SetFluidGasSolidViscosity.cpp
        FlowGravityForcing.h         
        LSLocateGasInterface.h       
        SetFluidGasSolidDensity.h    
        SetLSProperties.h
        GravityForcing.h             
        LevelSetInitialCondition.h   
        SetFluidGasSolidViscosity.h  
        TagLSRefinementCells.h)

FIND_PACKAGE(IBAMR 0.10.0 REQUIRED)
TARGET_LINK_LIBRARIES(main3d IBAMR::IBAMR3d)
# IBAMR saves the flags it used to compile - you can reuse them if you want to
SET(CMAKE_CXX_FLAGS ${IBAMR_CXX_FLAGS})
