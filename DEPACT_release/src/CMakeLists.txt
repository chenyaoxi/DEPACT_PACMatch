include_directories(${MYCODE_SOURCE_ROOT} ${OPENBABEL_BINARY_ROOT}/include/openbabel-2.0)
add_library(mycode samplerefdist.cpp atomtypessmarts.cpp fmmatches.cpp mmmatches.cpp extendssas.cpp bfdecoupling.cpp buildpocket.cpp statis_metal.cpp)
