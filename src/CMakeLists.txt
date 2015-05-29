
set(LIBBNET_SOURCE bnet.h bnet.c)
set(LIBBNET_INCLUDE bnet.h)

add_library(bnet SHARED ${LIBBNET_SOURCE})
target_link_libraries(bnet ${GSL_LIBRARIES})
install(TARGETS bnet DESTINATION lib)
install(FILES ${LIBBNET_INCLUDE} DESTINATION include)
