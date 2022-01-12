
SET ( SRC
    Input.cpp
    )

SET ( HEADERS
    Input.h
    Output.h
    VTK.h
    )

ADD_LIBRARY ( io  ${SRC} ${HEADERS} )



