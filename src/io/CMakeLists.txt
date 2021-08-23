
SET ( SRC
    Input.cpp
    )

SET ( HEADERS
    Input.h
    VTK.h
    vector_func.h
    Output.h
    )

ADD_LIBRARY ( io  ${SRC} ${HEADERS} )



