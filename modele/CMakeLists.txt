include_directories(${PROJECT_SOURCE_DIR}/slib)

find_package(TCLAP REQUIRED)
include_directories(${TCLAP_INCLUDE_DIR})

foreach (PRG
    giss2nc
    etopo1_ice overlap_hntr
    # make_topo oneway

    # Obsolete
    # make_topo_icebin icebin22m
)
    add_executable(${PRG} ${PRG}.cpp)
    target_link_libraries (${PRG} icebin ${EXTERNAL_LIBS})
    install(TARGETS ${PRG} DESTINATION bin)
endforeach()

add_executable(make_topo_f Z1QX1N.BS1.F)
target_link_libraries (make_topo_f icebin ${EXTERNAL_LIBS})
install(TARGETS make_topo_f DESTINATION bin)
