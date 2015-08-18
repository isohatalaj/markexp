
set(LIBBNET_SOURCE bnet.h bnet.c)
set(LIBBNET_INCLUDE bnet.h)

set(LIBMEXP_SOURCE mexp.h mexp.c)
set(LIBMEXP_INCLUDE mexp.h)

add_library(bnet SHARED ${LIBBNET_SOURCE})
target_link_libraries(bnet ${GSL_LIBRARIES})

add_library(mexp SHARED ${LIBMEXP_SOURCE})
target_link_libraries(mexp ${GSL_LIBRARIES})

install(TARGETS bnet DESTINATION lib)
install(FILES ${LIBBNET_INCLUDE} DESTINATION include)

install(TARGETS mexp DESTINATION lib)
install(FILES ${LIBMEXP_INCLUDE} DESTINATION include)
