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
CMAKE_SOURCE_DIR = /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build

# Include any dependencies generated for this target.
include test/CMakeFiles/DEPACTNativePocket.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/DEPACTNativePocket.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/DEPACTNativePocket.dir/flags.make

test/CMakeFiles/DEPACTNativePocket.dir/getnativepocket.o: test/CMakeFiles/DEPACTNativePocket.dir/flags.make
test/CMakeFiles/DEPACTNativePocket.dir/getnativepocket.o: ../test/getnativepocket.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/DEPACTNativePocket.dir/getnativepocket.o"
	cd /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DEPACTNativePocket.dir/getnativepocket.o -c /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/test/getnativepocket.cpp

test/CMakeFiles/DEPACTNativePocket.dir/getnativepocket.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DEPACTNativePocket.dir/getnativepocket.i"
	cd /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/test/getnativepocket.cpp > CMakeFiles/DEPACTNativePocket.dir/getnativepocket.i

test/CMakeFiles/DEPACTNativePocket.dir/getnativepocket.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DEPACTNativePocket.dir/getnativepocket.s"
	cd /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/test/getnativepocket.cpp -o CMakeFiles/DEPACTNativePocket.dir/getnativepocket.s

# Object files for target DEPACTNativePocket
DEPACTNativePocket_OBJECTS = \
"CMakeFiles/DEPACTNativePocket.dir/getnativepocket.o"

# External object files for target DEPACTNativePocket
DEPACTNativePocket_EXTERNAL_OBJECTS =

test/DEPACTNativePocket: test/CMakeFiles/DEPACTNativePocket.dir/getnativepocket.o
test/DEPACTNativePocket: test/CMakeFiles/DEPACTNativePocket.dir/build.make
test/DEPACTNativePocket: lib/libcombine.a
test/DEPACTNativePocket: test/CMakeFiles/DEPACTNativePocket.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable DEPACTNativePocket"
	cd /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DEPACTNativePocket.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/DEPACTNativePocket.dir/build: test/DEPACTNativePocket

.PHONY : test/CMakeFiles/DEPACTNativePocket.dir/build

test/CMakeFiles/DEPACTNativePocket.dir/clean:
	cd /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/test && $(CMAKE_COMMAND) -P CMakeFiles/DEPACTNativePocket.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/DEPACTNativePocket.dir/clean

test/CMakeFiles/DEPACTNativePocket.dir/depend:
	cd /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/test /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/test /home/cyx/workspace/hal/DEPACT_PACMatch/DEPACT_release/build/test/CMakeFiles/DEPACTNativePocket.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/DEPACTNativePocket.dir/depend

