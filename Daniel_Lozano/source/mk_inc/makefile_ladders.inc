# Include this to compile only the files needed for running the GS of a Ladder

# The directory where the object files will be created
OBJDIR = ../bin/obj

# The directory where the executable will be created
EXDIR := ../bin

# OBJMAIN gives the name of the main file to use for compiling
OBJMAIN := Ground_State_Ladder.o
DEPSMAIN := tntLadders.h

# SUBDIRS lists the subdirectories containing source code
SUBDIRS := functions

# List of files for each of the sub-directories
OBJstar := tntMpsCreateMpoLadder.o
DEPSstar := ../tntLadders.h

# Append the directory to the list of object files
OBJstar := $(patsubst %,functions/%,$(OBJstar))

# Add the files to the list of files to be compiled
OBJ := $(OBJMAIN) $(OBJstar)

# Append the object directory to the OBJ files.
OBJ := $(patsubst %,$(OBJDIR)/%,$(OBJ))
