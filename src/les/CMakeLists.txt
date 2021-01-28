SET ( MAINSRC

           main_LES.cpp
)


SET( executable   bdchmp )

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/lbsolver/")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/io/")
ADD_EXECUTABLE(        ${executable} ${MAINSRC} )
TARGET_LINK_LIBRARIES( ${executable} lbsolver io ${MPI_LIBRARIES})                                   