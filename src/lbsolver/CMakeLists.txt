
SET ( SRC
    Geo.cpp
    LBfield.cpp
    LBgeometry.cpp
    LBsnippets.cpp
    )

SET ( HEADERS
    Field.h
    Geo.h
    LBbndmpi.h
    LBbounceback.h
    LBboundary.h
    LBcollision.h
    LBcollision2phase.h
    LBd2q9.h
    LBd3q19.h
    LBfield.h
    LBfreeFlowCartesian.h
    LBfreeSlipCartesian.h
    LBfreeslipsolid.h
    LBgeometry.h
    LBglobal.h
    LBgrid.h
    LBhalfwaybb.h
    LBinitiatefield.h
    LBinletoutlet.h
    LBlatticetypes.h
    LBmacroscopic.h
    LBmonlatmpi.h
    LBnodes.h
    LBpressurebnd.h
    LBsnippets.h
    LButilities.h
    LBvtk.h
    LBrheology.h
    defines.h
    )

ADD_LIBRARY ( lbsolver  ${SRC} ${HEADERS} )



