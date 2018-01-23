add_executable(cardiac_ecc_fortran cardiac_ecc.F90)

target_link_libraries(cardiac_ecc_fortran PRIVATE opencmisslibs)

install(TARGETS cardiac_ecc_fortran DESTINATION bin)
