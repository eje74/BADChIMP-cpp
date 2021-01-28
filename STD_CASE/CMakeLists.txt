SET ( MAINSRC

           main_new_geo.cpp
)


SET( executable   bdchmp )

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/")
ADD_EXECUTABLE(        ${executable} ${MAINSRC} )
TARGET_LINK_LIBRARIES( ${executable} LBSOLVER ${MPI_LIBRARIES})                                   