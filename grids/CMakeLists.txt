include_directories(${PROJECT_SOURCE_DIR}/slib)


foreach (GRID
	modele_ll
	searise
	overlap

	# =========== Grid generators: Spherical Grids (for the GCM)
#	modele_ll_g2x2_5	# ModelE 2x2.5 grid, in vicinity of Greenland
#	modele_ll_g4x5		# ModelE 4x5 grid, in vicinity of Greenland
#	ga_2x2_5			# ModelE 2x2.5 grid, in vicinity of Greenland and Antarctica

	# =========== Grid generators: Cartesion Grids (for the Ice Model)
#	searise100			# 100km grid over ?<continent>? with SeaRISE projection
#	searise50			# 50km grid over ?<continent>? with SeaRISE projection
#	searise_g20			# 20km grid over Greenland with SeaRISE projection
#	searise_a			# The actual 5km grid used by SeaRISE over Antarctica
#	searise_g			# The actual 5km grid used by SeaRISE over Greenland
#	hires				# A 1km grid over Greenland with SeaRISE projection
)
	add_executable(${GRID} ${GRID}.cpp)
	target_link_libraries (${GRID} icebin ${EXTERNAL_LIBS})
    install(TARGETS ${GRID} DESTINATION bin)
endforeach()
