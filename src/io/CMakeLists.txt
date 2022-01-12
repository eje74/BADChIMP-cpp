
#SET ( SRC )

SET ( HEADERS
    Input.h
    VTK.h
    vector_func.h
    Output.h
    VTK.h
    )

#ADD_LIBRARY ( io  ${SRC} ${HEADERS} )
ADD_LIBRARY ( io  ${HEADERS} )
set_target_properties(io PROPERTIES LINKER_LANGUAGE CXX)
