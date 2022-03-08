
#SET ( SRC )

SET ( HEADERS
    Input.h
    Output.h
    VTK.h
    )

#ADD_LIBRARY ( io  ${SRC} ${HEADERS} )
#ADD_LIBRARY ( io  ${HEADERS} )
### This is a header-only target (no source-files) and must be included like this
add_library( io INTERFACE )
target_include_directories( io INTERFACE /. )
