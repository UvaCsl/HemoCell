# Add all `*.cpp` files to library `${TARGET}`. The libraries are marked as
# `EXLUDE_FROM_ALL`, such that they are only included in the make targets when
# when required.
function(add_source_to_library TARGET DIRLIST)
        foreach(d IN LISTS DIRLIST)
                file(GLOB_RECURSE SRC "${d}/*.cpp")
                list(APPEND SRC_FILES ${SRC})
        endforeach()
        add_library(${TARGET} STATIC EXCLUDE_FROM_ALL ${SRC_FILES})
endfunction(add_source_to_library)
