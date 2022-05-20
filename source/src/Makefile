prefix = ${GLB_DIR}
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

CC = g++
CFLAGS = -Wall
local_CFLAGS = -g -O4
debug_CFLAGS = -DDEBUG
shared_CFLAGS = -fPIC -shared -Wl,-soname,

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)



BIN = supernova 
OBJ = detector.o pinched.o pysnow.o main.o 
SONAME = libsnowflux.so.1.1.0
SOBJ = libsnowflux.so

all: $(BIN)



supernova: main.o detector.o pinched.o pysnow.o
	$(CC) main.o detector.o pinched.o pysnow.o -o supernova  $(LDFLAGS) $(local_LDFLAGS)

%.o : %.cpp
	$(CC) $(debug_CFLAGS) $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)

shared:
	$(CC) $(debug_CFLAGS) $(CFLAGS) $(local_CFLAGS) $(shared_CFLAGS)$(SONAME) -o $(SOBJ) *.cpp $(INCFLAGS) $(local_LDFLAGS)

install:
	cp $(BIN) ../

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ) $(SOBJ)
