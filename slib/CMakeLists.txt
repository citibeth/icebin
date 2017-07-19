include_directories(${PROJECT_SOURCE_DIR}/slib)

# To turn of GCC errors in one file...
# http://stackoverflow.com/questions/59670/how-to-get-rid-of-deprecated-conversion-from-string-constant-to-char-warnin/541079#541079
#pragma GCC diagnostic ignored "-Wwrite-strings"


#message("XXX0 " ${PROJECT_SOURCE_DIR}/slib)
#message("XXX1 " ${GMP_INCLUDE_DIR})
#message("XXX2 " ${CGAL_INCLUDE_DIR})
#message("XXX3 " ${BLITZ_INCLUDE_DIR})
#message("XXX4 " ${Boost_INCLUDE_DIRS})
#message("XXX7 " ${BOOST_ROOT})
#message("XXX8 " ${PETSC_INCLUDES})
#message("XXX9 " ${NETCDF_INCLUDES})
#
#message("PROJECT_SOURCE_DIR " ${PROJECT_SOURCE_DIR})

# -----------------------------------------------------------
#     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/f90blitz.f90


add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/f90blitz_f.f90
        ${CMAKE_CURRENT_BINARY_DIR}/f90blitz_f.h
    COMMAND
    python ${IBMISC_INCLUDE_DIR}/ibmisc/f90blitz_f.py icebin_f90blitz
        ${CMAKE_CURRENT_BINARY_DIR}/f90blitz_f.f90
        ${CMAKE_CURRENT_BINARY_DIR}/f90blitz_f.h)

#add_custom_target(makeme ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/f90blitz.f90)

set(icebin_SOURCES
    icebin/error.cpp
    icebin/Grid.cpp
    icebin/IceRegridder.cpp
    icebin/smoother.cpp
    icebin/GCMRegridder.cpp
    icebin/IceRegridder_L0.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/f90blitz_f.f90
)

if (BUILD_COUPLER)
    list(APPEND icebin_SOURCES
        # Coupler...
        icebin/multivec.cpp
        icebin/GCMCoupler.cpp
        icebin/IceCoupler.cpp
        icebin/contracts/contracts.cpp
    )

    if (BUILD_MODELE)
        list(APPEND icebin_SOURCES
            icebin/modele/api_f.f90
            icebin/modele/GCMCoupler_ModelE.cpp
            icebin/modele/hntr.cpp
            icebin/modele/z1qx1n_bs1.cpp
            icebin/modele/HNTR4.F90
        )

        if (BUILD_GRIDGEN)
            list(APPEND icebin_SOURCES
                icebin/modele/GridSpec_Hntr.cpp
                icebin/modele/clippers.cpp
            )
        endif()

    endif()

    if (USE_PISM)
        list(APPEND icebin_SOURCES
            icebin/pism/IceCoupler_PISM.cpp
        )
        if (BUILD_MODELE)
            list(APPEND icebin_SOURCES
                icebin/contracts/modele_pism.cpp
            )
        endif()
    endif()
endif()

if (BUILD_GRIDGEN)
    list(APPEND icebin_SOURCES
        icebin/gridgen/gridutil.cpp
        icebin/gridgen/clippers.cpp
        icebin/gridgen/GridSpec_LonLat.cpp
        icebin/gridgen/GridSpec_XY.cpp
        icebin/gridgen/GridSpec_Exchange.cpp
    )
endif()


if (USE_EVERYTRACE)
    if (BUILD_COUPLER)
      add_library (icebin SHARED ${EVERYTRACE_cf_mpi_REFADDR} ${icebin_SOURCES})
    else()
      add_library (icebin SHARED ${EVERYTRACE_cf_REFADDR} ${icebin_SOURCES})
    endif()
else()
    add_library (icebin SHARED ${icebin_SOURCES})
endif()

target_link_libraries(icebin ${EXTERNAL_LIBS})
if (USE_EVERYTRACE)
    target_link_libraries(icebin ${EVERYTRACE_LIBRARY})
endif()

install(TARGETS icebin DESTINATION lib)



# install Icebin headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/f90blitz_f.h
    DESTINATION include/icebin)

install(DIRECTORY ${PROJECT_BINARY_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.mod")
