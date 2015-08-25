
set(LIBMEXP_SOURCE mexp.h mexp.c)
set(LIBMEXP_INCLUDE mexp.h)

add_library(mexp SHARED ${LIBMEXP_SOURCE})
target_link_libraries(mexp ${GSL_LIBRARIES})

install(TARGETS mexp DESTINATION lib)
install(FILES ${LIBMEXP_INCLUDE} DESTINATION include)
