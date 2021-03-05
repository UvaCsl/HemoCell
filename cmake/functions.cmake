# add all `*.cpp` files to library `${TARGET}`
function(add_source_to_library TARGET DIRLIST)
        foreach(d IN LISTS DIRLIST)
                file(GLOB_RECURSE SRC "${d}/*.cpp")
                list(APPEND SRC_FILES ${SRC})
        endforeach()
        add_library(${TARGET} STATIC ${SRC_FILES})
endfunction(add_source_to_library)

# create executable from a directory (./example, ./cases)
# add source files with `*.cpp` or `*.h*` to executable `${TARGET}`
function(add_executable_from_source TARGET DIRLIST)
        foreach(d IN LISTS DIRLIST)
                file(GLOB_RECURSE SRC "${d}/*.cpp" "${d}/*.h*")
                list(APPEND SRC_FILES ${SRC})
        endforeach()
        add_executable(${TARGET} ${SRC_FILES})
endfunction(add_executable_from_source)
