include_directories(${PROJECT_SOURCE_DIR}/slib)

foreach (PRG
    giss2nc
)
    add_executable(${PRG} ${PRG}.cpp)
    target_link_libraries (${PRG} icebin ${EXTERNAL_LIBS})
    install(TARGETS ${PRG} DESTINATION bin)
endforeach()
