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
CMAKE_SOURCE_DIR = /home/cyx/workspace/DEPACT_release

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cyx/workspace/DEPACT_release/build

# Include any dependencies generated for this target.
include test/CMakeFiles/DEPACTPocketDesign.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/DEPACTPocketDesign.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/DEPACTPocketDesign.dir/flags.make

test/CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.o: test/CMakeFiles/DEPACTPocketDesign.dir/flags.make
test/CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.o: ../test/denovopocketdesign.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/DEPACT_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.o"
	cd /home/cyx/workspace/DEPACT_release/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.o -c /home/cyx/workspace/DEPACT_release/test/denovopocketdesign.cpp

test/CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.i"
	cd /home/cyx/workspace/DEPACT_release/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/DEPACT_release/test/denovopocketdesign.cpp > CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.i

test/CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.s"
	cd /home/cyx/workspace/DEPACT_release/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/DEPACT_release/test/denovopocketdesign.cpp -o CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.s

# Object files for target DEPACTPocketDesign
DEPACTPocketDesign_OBJECTS = \
"CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.o"

# External object files for target DEPACTPocketDesign
DEPACTPocketDesign_EXTERNAL_OBJECTS =

test/DEPACTPocketDesign: test/CMakeFiles/DEPACTPocketDesign.dir/denovopocketdesign.o
test/DEPACTPocketDesign: test/CMakeFiles/DEPACTPocketDesign.dir/build.make
test/DEPACTPocketDesign: lib/libcombine.a
test/DEPACTPocketDesign: test/CMakeFiles/DEPACTPocketDesign.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/cyx/workspace/DEPACT_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable DEPACTPocketDesign"
	cd /home/cyx/workspace/DEPACT_release/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DEPACTPocketDesign.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/DEPACTPocketDesign.dir/build: test/DEPACTPocketDesign

.PHONY : test/CMakeFiles/DEPACTPocketDesign.dir/build

test/CMakeFiles/DEPACTPocketDesign.dir/clean:
	cd /home/cyx/workspace/DEPACT_release/build/test && $(CMAKE_COMMAND) -P CMakeFiles/DEPACTPocketDesign.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/DEPACTPocketDesign.dir/clean

test/CMakeFiles/DEPACTPocketDesign.dir/depend:
	cd /home/cyx/workspace/DEPACT_release/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cyx/workspace/DEPACT_release /home/cyx/workspace/DEPACT_release/test /home/cyx/workspace/DEPACT_release/build /home/cyx/workspace/DEPACT_release/build/test /home/cyx/workspace/DEPACT_release/build/test/CMakeFiles/DEPACTPocketDesign.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/DEPACTPocketDesign.dir/depend
