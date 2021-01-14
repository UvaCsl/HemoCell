# add source code to directories
function(add_source_to_library TARGET DIRLIST)
        foreach(d IN LISTS DIRLIST)
                file(GLOB_RECURSE SRC "${d}/*.cpp")
                list(APPEND SRC_FILES ${SRC})
        endforeach()
        add_library(${TARGET} STATIC ${SRC_FILES})
endfunction(add_source_to_library)
