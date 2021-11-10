link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(testCS TestProteinRep.cpp)
add_executable(testScoring TestScoring.cpp)
add_executable(DesignSeq TestDesignSeq.cpp)
add_executable(design40 Test40PDB.cpp)
add_executable(design40Indv Test40PDBIndividual.cpp)

target_link_libraries (testCS loop)
target_link_libraries (testScoring loop)
target_link_libraries (DesignSeq loop)
target_link_libraries (design40 loop)
target_link_libraries (design40Indv loop)