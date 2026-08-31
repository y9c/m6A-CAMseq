# add_tag
# gcc -O3 -o add_tag add_tag.c -I/tools/include/ -L/tools/lib/ -lhts -lpthread -lz -lm -lbz2 -llzma
# Makefile for add_tag from add_tag.c

# Compiler
CC = gcc

# Executable name
TARGET = add_tag

# Source file
SRC = add_tag.c

# Compiler flags
# -O3 for optimization
# -Wall -Wextra for more warnings (good practice)
# -g for debugging information (can be removed for release builds if desired)
CFLAGS = -O3 -Wall -Wextra -g

# Include directories (adjust /tools/include if it's a variable or specific path)
# If HTSLIB_DIR is an environment variable or can be set:
# HTSLIB_DIR ?= /tools 
# INCLUDES = -I$(HTSLIB_DIR)/include
# For your specific case:
# INCLUDES = -I/tools/include/

# Library directories (adjust /tools/lib if it's a variable or specific path)
# LDFLAGS = -L$(HTSLIB_DIR)/lib
# For your specific case:
# LDFLAGS = -L/tools/lib/

# cd /home/yec/tools/htslib
# ./configure --prefix=/home/yec/tools/htslib_install 
# make
# make install
INCLUDES = -I/home/yec/tools/htslib_install/include
LDFLAGS = -L/home/yec/tools/htslib_install/lib

# Libraries to link against
# Order can sometimes matter for libraries
LIBS = -lhts -lpthread -lz -lm -lbz2 -llzma 

# Default target (what happens when you just type "make")
all: $(TARGET)

# Rule to link the executable
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SRC) $(LDFLAGS) $(LIBS)

# Rule to clean up compiled files
clean:
	rm -f $(TARGET) *.o core.* core

# Phony targets are targets that are not actual files
.PHONY: all clean
