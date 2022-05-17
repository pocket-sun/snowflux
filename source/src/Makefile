prefix = ${GLB_DIR}
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

CC = g++
CFLAGS = -Wall
local_CFLAGS = -O4
debug_CFLAGS = -DDEBUG

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)



BIN = supernova 
OBJ = main.o detector.o pinched.o

all: $(BIN)



supernova: main.o detector.o pinched.o
	$(CC) main.o detector.o pinched.o -o supernova  $(LDFLAGS) $(local_LDFLAGS)

%.o : %.cpp
	$(CC) $(debug_CFLAGS) $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)

install:
	cp $(BIN) ../

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)
