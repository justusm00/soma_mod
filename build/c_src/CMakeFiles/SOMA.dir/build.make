# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.24.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.24.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/justusmulthaup/soma_mod

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/justusmulthaup/soma_mod/build

# Include any dependencies generated for this target.
include c_src/CMakeFiles/SOMA.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include c_src/CMakeFiles/SOMA.dir/compiler_depend.make

# Include the progress variables for this target.
include c_src/CMakeFiles/SOMA.dir/progress.make

# Include the compile flags for this target's objects.
include c_src/CMakeFiles/SOMA.dir/flags.make

c_src/CMakeFiles/SOMA.dir/soma.c.o: c_src/CMakeFiles/SOMA.dir/flags.make
c_src/CMakeFiles/SOMA.dir/soma.c.o: /Users/justusmulthaup/soma_mod/c_src/soma.c
c_src/CMakeFiles/SOMA.dir/soma.c.o: c_src/CMakeFiles/SOMA.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/justusmulthaup/soma_mod/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object c_src/CMakeFiles/SOMA.dir/soma.c.o"
	cd /Users/justusmulthaup/soma_mod/build/c_src && /opt/homebrew/Cellar/llvm/15.0.3/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT c_src/CMakeFiles/SOMA.dir/soma.c.o -MF CMakeFiles/SOMA.dir/soma.c.o.d -o CMakeFiles/SOMA.dir/soma.c.o -c /Users/justusmulthaup/soma_mod/c_src/soma.c

c_src/CMakeFiles/SOMA.dir/soma.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SOMA.dir/soma.c.i"
	cd /Users/justusmulthaup/soma_mod/build/c_src && /opt/homebrew/Cellar/llvm/15.0.3/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/justusmulthaup/soma_mod/c_src/soma.c > CMakeFiles/SOMA.dir/soma.c.i

c_src/CMakeFiles/SOMA.dir/soma.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SOMA.dir/soma.c.s"
	cd /Users/justusmulthaup/soma_mod/build/c_src && /opt/homebrew/Cellar/llvm/15.0.3/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/justusmulthaup/soma_mod/c_src/soma.c -o CMakeFiles/SOMA.dir/soma.c.s

# Object files for target SOMA
SOMA_OBJECTS = \
"CMakeFiles/SOMA.dir/soma.c.o"

# External object files for target SOMA
SOMA_EXTERNAL_OBJECTS =

c_src/SOMA: c_src/CMakeFiles/SOMA.dir/soma.c.o
c_src/SOMA: c_src/CMakeFiles/SOMA.dir/build.make
c_src/SOMA: c_src/libsoma_lib.a
c_src/SOMA: /opt/homebrew/Cellar/open-mpi/4.1.4_2/lib/libmpi.dylib
c_src/SOMA: /opt/homebrew/Cellar/hdf5-mpi/1.12.2_1/lib/libhdf5.dylib
c_src/SOMA: /opt/homebrew/opt/libaec/lib/libsz.dylib
c_src/SOMA: /Library/Developer/CommandLineTools/SDKs/MacOSX13.0.sdk/usr/lib/libz.tbd
c_src/SOMA: /Library/Developer/CommandLineTools/SDKs/MacOSX13.0.sdk/usr/lib/libdl.tbd
c_src/SOMA: /Library/Developer/CommandLineTools/SDKs/MacOSX13.0.sdk/usr/lib/libm.tbd
c_src/SOMA: /opt/homebrew/Cellar/hdf5-mpi/1.12.2_1/lib/libhdf5.dylib
c_src/SOMA: /opt/homebrew/opt/libaec/lib/libsz.dylib
c_src/SOMA: /Library/Developer/CommandLineTools/SDKs/MacOSX13.0.sdk/usr/lib/libz.tbd
c_src/SOMA: /Library/Developer/CommandLineTools/SDKs/MacOSX13.0.sdk/usr/lib/libdl.tbd
c_src/SOMA: /Library/Developer/CommandLineTools/SDKs/MacOSX13.0.sdk/usr/lib/libm.tbd
c_src/SOMA: c_src/CMakeFiles/SOMA.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/justusmulthaup/soma_mod/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable SOMA"
	cd /Users/justusmulthaup/soma_mod/build/c_src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SOMA.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
c_src/CMakeFiles/SOMA.dir/build: c_src/SOMA
.PHONY : c_src/CMakeFiles/SOMA.dir/build

c_src/CMakeFiles/SOMA.dir/clean:
	cd /Users/justusmulthaup/soma_mod/build/c_src && $(CMAKE_COMMAND) -P CMakeFiles/SOMA.dir/cmake_clean.cmake
.PHONY : c_src/CMakeFiles/SOMA.dir/clean

c_src/CMakeFiles/SOMA.dir/depend:
	cd /Users/justusmulthaup/soma_mod/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/justusmulthaup/soma_mod /Users/justusmulthaup/soma_mod/c_src /Users/justusmulthaup/soma_mod/build /Users/justusmulthaup/soma_mod/build/c_src /Users/justusmulthaup/soma_mod/build/c_src/CMakeFiles/SOMA.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : c_src/CMakeFiles/SOMA.dir/depend

