
# Make sure the linker can find the Glint2 library once it is built.
link_directories (${PROJECT_BINARY_DIR}/slib)

# ============= General Utilites
# --------- Find the overlap of two grids
add_executable (overlap overlap.cpp)
target_link_libraries (overlap glint2 ${EXTERNAL_LIBS}) 

# --------- Replay recorded output of GCM into Glint2.  For one-way coupling.
add_executable (desm desm.cpp)
target_link_libraries (desm glint2 ${EXTERNAL_LIBS}) 

# =========== Grid generators: Spherical Grids (for the GCM)

# ModelE 2x2.5 grid, in vicinity of Greenland
add_executable(modele_ll_g2x2_5 modele_ll_g2x2_5.cpp)
target_link_libraries (modele_ll_g2x2_5 glint2 ${EXTERNAL_LIBS})

# ModelE 4x5 grid, in vicinity of Greenland
add_executable(modele_ll_g4x5 modele_ll_g4x5.cpp)
target_link_libraries (modele_ll_g4x5 glint2 ${EXTERNAL_LIBS})

# ModelE 2x2.5 grid, in vicinity of Greenland and Antarctica
add_executable(ga_2x2_5 ga_2x2_5.cpp)
target_link_libraries (ga_2x2_5 glint2 ${EXTERNAL_LIBS})

# =========== Grid generators: Cartesion Grids (for the Ice Model)

# 100km grid over ?<continent>? with SeaRISE projection
add_executable(searise100 searise100.cpp)
target_link_libraries (searise100 glint2 ${EXTERNAL_LIBS})

# 50km grid over ?<continent>? with SeaRISE projection
add_executable(searise50 searise50.cpp)
target_link_libraries (searise50 glint2 ${EXTERNAL_LIBS})

# 20km grid over Greenland with SeaRISE projection
add_executable(searise_g20 searise_g20.cpp)
target_link_libraries (searise_g20 glint2 ${EXTERNAL_LIBS})

# The actual 5km grid used by SeaRISE over Antarctica
add_executable(searise_a searise_a.cpp)
target_link_libraries (searise_a glint2 ${EXTERNAL_LIBS})

# The actual 5km grid used by SeaRISE over Greenland
add_executable(searise_g searise_g.cpp)
target_link_libraries (searise_g glint2 ${EXTERNAL_LIBS})

# A 1km grid over Greenland with SeaRISE projection
add_executable(hires hires.cpp)
target_link_libraries (hires glint2 ${EXTERNAL_LIBS})

# ============= Examples and Tests
# ---------- Test the API used by ModelE

# Doesn't currently compile.
#add_executable(apitest apitest.cpp)
#target_link_libraries (apitest glint2 ${EXTERNAL_LIBS})

# ---------- Some example with Blitz++
#add_executable(blitztest blitztest.cpp)
#target_link_libraries (blitztest glint2 ${EXTERNAL_LIBS})

# --------- The PISM main program, but compiled from with Glint2.
# If this works, it shows we've set up the build properly.
#add_executable(pismr pismr.cpp)        # Doesn't currently compile
#target_link_libraries (pismr glint2 ${EXTERNAL_LIBS})

## Unit test for sparse matrix multiplication
#add_executable(smulttest smulttest.cpp)
#target_link_libraries (smulttest glint2 ${EXTERNAL_LIBS})
#
#add_executable(test_grid1 test_grid1.cpp)
#target_link_libraries (test_grid1 glint2 ${EXTERNAL_LIBS})
#
#add_executable(test_grid2 test_grid2.cpp)
#target_link_libraries (test_grid2 glint2 ${EXTERNAL_LIBS})

# ================================================
# ================================================

install(TARGETS
#       apitest
        desm
        ga_2x2_5
        modele_ll_g2x2_5
        modele_ll_g4x5
        hires
        overlap
#       pismr
        searise_a
        searise_g
    DESTINATION bin)

