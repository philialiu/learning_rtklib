# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /home/philia/Applications/cmake-3.14.4-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/philia/Applications/cmake-3.14.4-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/philia/Documents/learning_rtklib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/philia/Documents/learning_rtklib/build

# Include any dependencies generated for this target.
include CMakeFiles/satpos_example.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/satpos_example.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/satpos_example.dir/flags.make

CMakeFiles/satpos_example.dir/satpos_example.cpp.o: CMakeFiles/satpos_example.dir/flags.make
CMakeFiles/satpos_example.dir/satpos_example.cpp.o: ../satpos_example.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/philia/Documents/learning_rtklib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/satpos_example.dir/satpos_example.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/satpos_example.dir/satpos_example.cpp.o -c /home/philia/Documents/learning_rtklib/satpos_example.cpp

CMakeFiles/satpos_example.dir/satpos_example.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/satpos_example.dir/satpos_example.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/philia/Documents/learning_rtklib/satpos_example.cpp > CMakeFiles/satpos_example.dir/satpos_example.cpp.i

CMakeFiles/satpos_example.dir/satpos_example.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/satpos_example.dir/satpos_example.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/philia/Documents/learning_rtklib/satpos_example.cpp -o CMakeFiles/satpos_example.dir/satpos_example.cpp.s

# Object files for target satpos_example
satpos_example_OBJECTS = \
"CMakeFiles/satpos_example.dir/satpos_example.cpp.o"

# External object files for target satpos_example
satpos_example_EXTERNAL_OBJECTS =

satpos_example: CMakeFiles/satpos_example.dir/satpos_example.cpp.o
satpos_example: CMakeFiles/satpos_example.dir/build.make
satpos_example: CMakeFiles/satpos_example.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/philia/Documents/learning_rtklib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable satpos_example"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/satpos_example.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/satpos_example.dir/build: satpos_example

.PHONY : CMakeFiles/satpos_example.dir/build

CMakeFiles/satpos_example.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/satpos_example.dir/cmake_clean.cmake
.PHONY : CMakeFiles/satpos_example.dir/clean

CMakeFiles/satpos_example.dir/depend:
	cd /home/philia/Documents/learning_rtklib/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/philia/Documents/learning_rtklib /home/philia/Documents/learning_rtklib /home/philia/Documents/learning_rtklib/build /home/philia/Documents/learning_rtklib/build /home/philia/Documents/learning_rtklib/build/CMakeFiles/satpos_example.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/satpos_example.dir/depend

