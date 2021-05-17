# Attaches the discovered PARMETIS directories to each target. Note: this is
# only done when `PARMETIS_FOUND` is set to true. Otherwise, the necessary
# libraries were not found and the targets do not require to be linked.
function(ConfigureParmetis TARGETS)
        if(${PARMETIS_FOUND})
                foreach(TARGET ${TARGETS})
                        target_include_directories(${TARGET} PRIVATE ${PARMETIS_INCLUDE_DIR})
                        target_include_directories(${TARGET} PRIVATE ${METIS_INCLUDE_DIR})
                        target_link_libraries(${TARGET} PRIVATE ${PARMETIS_LIBRARY})
                        target_link_libraries(${TARGET} PRIVATE ${METIS_LIBRARY})
                endforeach()
        endif()
endfunction(ConfigureParmetis)
