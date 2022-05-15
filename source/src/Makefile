prefix = ${GLB_DIR}
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

CC = g++
CFLAGS = -Wall
local_CFLAGS = -O4

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)



BIN = supernova 
OBJ = supernova.o pinched.o

all: $(BIN)



supernova: supernova.o pinched.o
	$(CC) supernova.o pinched.o -o supernova  $(LDFLAGS) $(local_LDFLAGS)

%.o : %.cpp
	$(CC) $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
	$(CC) $(CFLAGS) pinched.cpp -c

install:
	cp $(BIN) ../

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)
