# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cyx/workspace/PACMatch_release

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cyx/workspace/PACMatch_release/build

# Include any dependencies generated for this target.
include proteinrep/src/CMakeFiles/nspproteinrep.dir/depend.make

# Include the progress variables for this target.
include proteinrep/src/CMakeFiles/nspproteinrep.dir/progress.make

# Include the compile flags for this target's objects.
include proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make

proteinrep/src/CMakeFiles/nspproteinrep.dir/intatomkey.o: proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make
proteinrep/src/CMakeFiles/nspproteinrep.dir/intatomkey.o: ../proteinrep/src/intatomkey.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object proteinrep/src/CMakeFiles/nspproteinrep.dir/intatomkey.o"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nspproteinrep.dir/intatomkey.o -c /home/cyx/workspace/PACMatch_release/proteinrep/src/intatomkey.cpp

proteinrep/src/CMakeFiles/nspproteinrep.dir/intatomkey.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nspproteinrep.dir/intatomkey.i"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/PACMatch_release/proteinrep/src/intatomkey.cpp > CMakeFiles/nspproteinrep.dir/intatomkey.i

proteinrep/src/CMakeFiles/nspproteinrep.dir/intatomkey.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nspproteinrep.dir/intatomkey.s"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/PACMatch_release/proteinrep/src/intatomkey.cpp -o CMakeFiles/nspproteinrep.dir/intatomkey.s

proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbrecord.o: proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make
proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbrecord.o: ../proteinrep/src/pdbrecord.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbrecord.o"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nspproteinrep.dir/pdbrecord.o -c /home/cyx/workspace/PACMatch_release/proteinrep/src/pdbrecord.cpp

proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbrecord.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nspproteinrep.dir/pdbrecord.i"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/PACMatch_release/proteinrep/src/pdbrecord.cpp > CMakeFiles/nspproteinrep.dir/pdbrecord.i

proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbrecord.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nspproteinrep.dir/pdbrecord.s"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/PACMatch_release/proteinrep/src/pdbrecord.cpp -o CMakeFiles/nspproteinrep.dir/pdbrecord.s

proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbreader.o: proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make
proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbreader.o: ../proteinrep/src/pdbreader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbreader.o"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nspproteinrep.dir/pdbreader.o -c /home/cyx/workspace/PACMatch_release/proteinrep/src/pdbreader.cpp

proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbreader.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nspproteinrep.dir/pdbreader.i"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/PACMatch_release/proteinrep/src/pdbreader.cpp > CMakeFiles/nspproteinrep.dir/pdbreader.i

proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbreader.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nspproteinrep.dir/pdbreader.s"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/PACMatch_release/proteinrep/src/pdbreader.cpp -o CMakeFiles/nspproteinrep.dir/pdbreader.s

proteinrep/src/CMakeFiles/nspproteinrep.dir/aminoacidseq.o: proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make
proteinrep/src/CMakeFiles/nspproteinrep.dir/aminoacidseq.o: ../proteinrep/src/aminoacidseq.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object proteinrep/src/CMakeFiles/nspproteinrep.dir/aminoacidseq.o"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nspproteinrep.dir/aminoacidseq.o -c /home/cyx/workspace/PACMatch_release/proteinrep/src/aminoacidseq.cpp

proteinrep/src/CMakeFiles/nspproteinrep.dir/aminoacidseq.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nspproteinrep.dir/aminoacidseq.i"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/PACMatch_release/proteinrep/src/aminoacidseq.cpp > CMakeFiles/nspproteinrep.dir/aminoacidseq.i

proteinrep/src/CMakeFiles/nspproteinrep.dir/aminoacidseq.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nspproteinrep.dir/aminoacidseq.s"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/PACMatch_release/proteinrep/src/aminoacidseq.cpp -o CMakeFiles/nspproteinrep.dir/aminoacidseq.s

proteinrep/src/CMakeFiles/nspproteinrep.dir/aaconformer.o: proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make
proteinrep/src/CMakeFiles/nspproteinrep.dir/aaconformer.o: ../proteinrep/src/aaconformer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object proteinrep/src/CMakeFiles/nspproteinrep.dir/aaconformer.o"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nspproteinrep.dir/aaconformer.o -c /home/cyx/workspace/PACMatch_release/proteinrep/src/aaconformer.cpp

proteinrep/src/CMakeFiles/nspproteinrep.dir/aaconformer.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nspproteinrep.dir/aaconformer.i"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/PACMatch_release/proteinrep/src/aaconformer.cpp > CMakeFiles/nspproteinrep.dir/aaconformer.i

proteinrep/src/CMakeFiles/nspproteinrep.dir/aaconformer.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nspproteinrep.dir/aaconformer.s"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/PACMatch_release/proteinrep/src/aaconformer.cpp -o CMakeFiles/nspproteinrep.dir/aaconformer.s

proteinrep/src/CMakeFiles/nspproteinrep.dir/idealgeometries.o: proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make
proteinrep/src/CMakeFiles/nspproteinrep.dir/idealgeometries.o: ../proteinrep/src/idealgeometries.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object proteinrep/src/CMakeFiles/nspproteinrep.dir/idealgeometries.o"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nspproteinrep.dir/idealgeometries.o -c /home/cyx/workspace/PACMatch_release/proteinrep/src/idealgeometries.cpp

proteinrep/src/CMakeFiles/nspproteinrep.dir/idealgeometries.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nspproteinrep.dir/idealgeometries.i"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/PACMatch_release/proteinrep/src/idealgeometries.cpp > CMakeFiles/nspproteinrep.dir/idealgeometries.i

proteinrep/src/CMakeFiles/nspproteinrep.dir/idealgeometries.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nspproteinrep.dir/idealgeometries.s"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/PACMatch_release/proteinrep/src/idealgeometries.cpp -o CMakeFiles/nspproteinrep.dir/idealgeometries.s

proteinrep/src/CMakeFiles/nspproteinrep.dir/chaintree.o: proteinrep/src/CMakeFiles/nspproteinrep.dir/flags.make
proteinrep/src/CMakeFiles/nspproteinrep.dir/chaintree.o: ../proteinrep/src/chaintree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object proteinrep/src/CMakeFiles/nspproteinrep.dir/chaintree.o"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nspproteinrep.dir/chaintree.o -c /home/cyx/workspace/PACMatch_release/proteinrep/src/chaintree.cpp

proteinrep/src/CMakeFiles/nspproteinrep.dir/chaintree.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nspproteinrep.dir/chaintree.i"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/PACMatch_release/proteinrep/src/chaintree.cpp > CMakeFiles/nspproteinrep.dir/chaintree.i

proteinrep/src/CMakeFiles/nspproteinrep.dir/chaintree.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nspproteinrep.dir/chaintree.s"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/PACMatch_release/proteinrep/src/chaintree.cpp -o CMakeFiles/nspproteinrep.dir/chaintree.s

# Object files for target nspproteinrep
nspproteinrep_OBJECTS = \
"CMakeFiles/nspproteinrep.dir/intatomkey.o" \
"CMakeFiles/nspproteinrep.dir/pdbrecord.o" \
"CMakeFiles/nspproteinrep.dir/pdbreader.o" \
"CMakeFiles/nspproteinrep.dir/aminoacidseq.o" \
"CMakeFiles/nspproteinrep.dir/aaconformer.o" \
"CMakeFiles/nspproteinrep.dir/idealgeometries.o" \
"CMakeFiles/nspproteinrep.dir/chaintree.o"

# External object files for target nspproteinrep
nspproteinrep_EXTERNAL_OBJECTS =

lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/intatomkey.o
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbrecord.o
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/pdbreader.o
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/aminoacidseq.o
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/aaconformer.o
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/idealgeometries.o
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/chaintree.o
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/build.make
lib/libnspproteinrep.a: proteinrep/src/CMakeFiles/nspproteinrep.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/cyx/workspace/PACMatch_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX static library ../../lib/libnspproteinrep.a"
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && $(CMAKE_COMMAND) -P CMakeFiles/nspproteinrep.dir/cmake_clean_target.cmake
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nspproteinrep.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
proteinrep/src/CMakeFiles/nspproteinrep.dir/build: lib/libnspproteinrep.a

.PHONY : proteinrep/src/CMakeFiles/nspproteinrep.dir/build

proteinrep/src/CMakeFiles/nspproteinrep.dir/clean:
	cd /home/cyx/workspace/PACMatch_release/build/proteinrep/src && $(CMAKE_COMMAND) -P CMakeFiles/nspproteinrep.dir/cmake_clean.cmake
.PHONY : proteinrep/src/CMakeFiles/nspproteinrep.dir/clean

proteinrep/src/CMakeFiles/nspproteinrep.dir/depend:
	cd /home/cyx/workspace/PACMatch_release/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cyx/workspace/PACMatch_release /home/cyx/workspace/PACMatch_release/proteinrep/src /home/cyx/workspace/PACMatch_release/build /home/cyx/workspace/PACMatch_release/build/proteinrep/src /home/cyx/workspace/PACMatch_release/build/proteinrep/src/CMakeFiles/nspproteinrep.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : proteinrep/src/CMakeFiles/nspproteinrep.dir/depend
