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


# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.

set(icebin_SOURCES
    icebin/error.cpp
    icebin/Grid.cpp
    icebin/GCMRegridder.cpp
    icebin/IceRegridder_L0.cpp
)

if (BUILD_COUPLER)
    list(APPEND icebin_SOURCES
        # Coupler...
        icebin/GCMParams.cpp
        icebin/GCMCoupler.cpp
        icebin/IceModel_Writer.cpp
        icebin/IceModel_Decode.cpp
#       icebin/IceModel_DISMAL.cpp
        icebin/contracts/contracts.cpp
    )
endif()

if (USE_PISM)
    list(APPEND icebin_SOURCES
        icebin/pism/IceModel_PISM.cpp
        icebin/contracts/modele_pism.cpp
    )
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


# set(icebin_SOURCES 
#   libicebin_refaddr.cpp
#   galahad/qpt_c.cpp
#   galahad/zd11_c.cpp
#   giss/DynamicEnum.cpp
#   giss/IndexTranslator.cpp
#   giss/IndexTranslator2.cpp
#   giss/Proj2.cpp
#   giss/SparseMatrix.cpp
#   giss/cfnames.cpp
#   giss/cfnames_data.cpp
#   giss/geodesy.cpp
#   giss/ncutil.cpp
#   giss/sparsemult.cpp
#   giss/udunits2.cpp
#   giss/exit.cpp
#   giss/exit_f.F90
#   giss/VarTransformer.cpp
#   icebin/contracts/contracts.cpp
#   icebin/contracts/modele_pism.cpp
#   icebin/contracts/modele_dismal.cpp
#   giss/ConstantSet.cpp
#   giss/CouplingContract.cpp
#   icebin/ExchangeGrid.cpp
#   icebin/GCMCoupler.cpp
#   icebin/GCMParams.cpp
#   icebin/Grid.cpp
#   icebin/GridDomain.cpp
#   icebin/GridDomain_Listed.cpp
#   icebin/Grid_LonLat.cpp
#   icebin/Grid_XY.cpp
#   icebin/HCIndex.cpp
#   icebin/HeightClassifier.cpp
#   icebin/IceModel.cpp
#   icebin/IceModel_DISMAL.cpp
#   icebin/IceModel_Writer.cpp
#   icebin/IceModel_Decode.cpp
#   icebin/IceSheet.cpp
#   icebin/IceSheet_L0.cpp
#   icebin/MatrixMaker.cpp
#   icebin/MultiMatrix.cpp
#   icebin/clippers.cpp
#   icebin/gridutil.cpp
#   icebin/matrix_ops.cpp
#   icebin/contracts/modele_pism.cpp
#   icebin/modele/icebin_modele.cpp
#   icebin/modele/GCMCoupler_ModelE.cpp
#   icebin/modele/grids_ll.cpp
#   icebin/pism/VecBundleWriter.cpp
#   icebin/pism/IceModel_PISM.cpp
#   icebin/pism/NullTransportHydrology.cpp
#   icebin/pism/MassEnergyBudget.cpp
#   icebin/pism/PISMIceModel.cpp
#   icebin/pism/PSConstantICEBIN.cpp
#   icebin/pism/ICEBINEnthalpyConverter.cpp
#   icebin/read_grid.cpp
#   icebin/read_icemodel.cpp
#   icebin/read_matrixmaker.cpp
#   giss/f90blitz_f.f90
#   icebin/modele/icebin_modele_f.f90
#   # GALAHAD sources: turn off if USE_GALAHAD
#   galahad/eqp_f.f90
#   galahad/ncutil_f.f90
#   galahad/qpt_f.f90
#   galahad/x_loc_f.f90
#   galahad/zd11_f.f90
#   galahad/zd11_c.cpp
#   galahad/qpt_c.cpp)

add_library (icebin SHARED ${icebin_SOURCES} ${EVERYTRACE_REFADDR})

target_link_libraries(icebin ${EXTERNAL_LIBS})

install(TARGETS icebin DESTINATION lib)

# install Icebin headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")

install(DIRECTORY ${PROJECT_BINARY_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.mod")
