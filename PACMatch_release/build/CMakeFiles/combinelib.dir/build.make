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

# Utility rule file for combinelib.

# Include the progress variables for this target.
include CMakeFiles/combinelib.dir/progress.make

CMakeFiles/combinelib: ../backbone
CMakeFiles/combinelib: lib/libnspdataio.a
CMakeFiles/combinelib: lib/libnspdstl.a
CMakeFiles/combinelib: lib/libnspgeometry.a
CMakeFiles/combinelib: lib/libnspproteinrep.a
CMakeFiles/combinelib: lib/libnspdesignseq.a
CMakeFiles/combinelib: lib/libnoob.a
	cd /home/cyx/workspace/PACMatch_release/build/lib && rm -f libcombine.a
	cd /home/cyx/workspace/PACMatch_release/build/lib && ar rcT libcombine.a libnspdataio.a libnspdstl.a libnspgeometry.a libnspproteinrep.a libnspdesignseq.a libnoob.a
	cd /home/cyx/workspace/PACMatch_release/build/lib && ranlib libcombine.a

combinelib: CMakeFiles/combinelib
combinelib: CMakeFiles/combinelib.dir/build.make

.PHONY : combinelib

# Rule to build all files generated by this target.
CMakeFiles/combinelib.dir/build: combinelib

.PHONY : CMakeFiles/combinelib.dir/build

CMakeFiles/combinelib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/combinelib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/combinelib.dir/clean

CMakeFiles/combinelib.dir/depend:
	cd /home/cyx/workspace/PACMatch_release/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cyx/workspace/PACMatch_release /home/cyx/workspace/PACMatch_release /home/cyx/workspace/PACMatch_release/build /home/cyx/workspace/PACMatch_release/build /home/cyx/workspace/PACMatch_release/build/CMakeFiles/combinelib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/combinelib.dir/depend

