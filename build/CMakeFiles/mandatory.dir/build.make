# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/panda1/numerical-methods

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/panda1/numerical-methods/build

# Include any dependencies generated for this target.
include CMakeFiles/mandatory.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mandatory.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mandatory.dir/flags.make

CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o: CMakeFiles/mandatory.dir/flags.make
CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o: ../src/mandatory/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/panda1/numerical-methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o -c /home/panda1/numerical-methods/src/mandatory/main.cpp

CMakeFiles/mandatory.dir/src/mandatory/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mandatory.dir/src/mandatory/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/panda1/numerical-methods/src/mandatory/main.cpp > CMakeFiles/mandatory.dir/src/mandatory/main.cpp.i

CMakeFiles/mandatory.dir/src/mandatory/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mandatory.dir/src/mandatory/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/panda1/numerical-methods/src/mandatory/main.cpp -o CMakeFiles/mandatory.dir/src/mandatory/main.cpp.s

CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.requires:

.PHONY : CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.requires

CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.provides: CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/mandatory.dir/build.make CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.provides.build
.PHONY : CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.provides

CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.provides.build: CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o


# Object files for target mandatory
mandatory_OBJECTS = \
"CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o"

# External object files for target mandatory
mandatory_EXTERNAL_OBJECTS =

mandatory: CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o
mandatory: CMakeFiles/mandatory.dir/build.make
mandatory: CMakeFiles/mandatory.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/panda1/numerical-methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mandatory"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mandatory.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mandatory.dir/build: mandatory

.PHONY : CMakeFiles/mandatory.dir/build

CMakeFiles/mandatory.dir/requires: CMakeFiles/mandatory.dir/src/mandatory/main.cpp.o.requires

.PHONY : CMakeFiles/mandatory.dir/requires

CMakeFiles/mandatory.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mandatory.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mandatory.dir/clean

CMakeFiles/mandatory.dir/depend:
	cd /home/panda1/numerical-methods/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/panda1/numerical-methods /home/panda1/numerical-methods /home/panda1/numerical-methods/build /home/panda1/numerical-methods/build /home/panda1/numerical-methods/build/CMakeFiles/mandatory.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mandatory.dir/depend

