# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/lehatrutenb/Desktop/prog/fluid1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/lehatrutenb/Desktop/prog/fluid1/fluid_build

# Include any dependencies generated for this target.
include CMakeFiles/Fluid3Project.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Fluid3Project.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Fluid3Project.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Fluid3Project.dir/flags.make

CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o: CMakeFiles/Fluid3Project.dir/flags.make
CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o: /Users/lehatrutenb/Desktop/prog/fluid1/fluid3/fluid.cpp
CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o: CMakeFiles/Fluid3Project.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/lehatrutenb/Desktop/prog/fluid1/fluid_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o -MF CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o.d -o CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o -c /Users/lehatrutenb/Desktop/prog/fluid1/fluid3/fluid.cpp

CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/lehatrutenb/Desktop/prog/fluid1/fluid3/fluid.cpp > CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.i

CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/lehatrutenb/Desktop/prog/fluid1/fluid3/fluid.cpp -o CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.s

# Object files for target Fluid3Project
Fluid3Project_OBJECTS = \
"CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o"

# External object files for target Fluid3Project
Fluid3Project_EXTERNAL_OBJECTS =

Fluid3Project: CMakeFiles/Fluid3Project.dir/fluid3/fluid.cpp.o
Fluid3Project: CMakeFiles/Fluid3Project.dir/build.make
Fluid3Project: CMakeFiles/Fluid3Project.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/lehatrutenb/Desktop/prog/fluid1/fluid_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Fluid3Project"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Fluid3Project.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Fluid3Project.dir/build: Fluid3Project
.PHONY : CMakeFiles/Fluid3Project.dir/build

CMakeFiles/Fluid3Project.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Fluid3Project.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Fluid3Project.dir/clean

CMakeFiles/Fluid3Project.dir/depend:
	cd /Users/lehatrutenb/Desktop/prog/fluid1/fluid_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/lehatrutenb/Desktop/prog/fluid1 /Users/lehatrutenb/Desktop/prog/fluid1 /Users/lehatrutenb/Desktop/prog/fluid1/fluid_build /Users/lehatrutenb/Desktop/prog/fluid1/fluid_build /Users/lehatrutenb/Desktop/prog/fluid1/fluid_build/CMakeFiles/Fluid3Project.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/Fluid3Project.dir/depend

