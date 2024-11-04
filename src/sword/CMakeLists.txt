set(SWORD_SOURCES
    chain.cpp
    chain.hpp
    evalue.cpp
    evalue.hpp
    hash.cpp
    hash.hpp
    kmers.cpp
    kmers.hpp
    reader.cpp
    reader.hpp
    score_matrix.cpp
    score_matrix.hpp
)




foreach(_file IN ITEMS ${SWORD_SOURCES})
    if(EXISTS ${PROJECT_SOURCE_DIR}/patches/${_file}.patch)
        add_custom_command(
            OUTPUT 
                ${_file}
            COMMENT
                "Patching ${_file}"
            COMMAND
                ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py 
                    --input ${PROJECT_SOURCE_DIR}/vendor/sword/src/${_file} 
                    --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                    --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/sword/src/${_file}
                ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
        )
    else()
        add_custom_command(
            OUTPUT 
                ${_file}
            COMMAND
                cmake -E copy ${PROJECT_SOURCE_DIR}/vendor/sword/src/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/sword/src/${_file}
        )
    endif()
    set(SWORD_PATCHED_SOURCES ${SWORD_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

add_library(sword STATIC ${SWORD_PATCHED_SOURCES})
target_include_directories(sword PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(sword PUBLIC ${PROJECT_SOURCE_DIR}/vendor/sword/vendor/opal/src)