UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CXX = g++
	GDB = gdb
	LIBS = -lglui -lglut -lGLEW -lGL -lGLU -larmadillo -lCGAL -lgmp -frounding-math
	BADLINKS =
	OPTFLAGS = -w
endif
ifeq ($(UNAME_S),Darwin)
	CXX = clang++
	GDB = lldb
	LIBS = -framework GLUI -framework OpenGL -framework GLUT -lGLEW -larmadillo
	BADLINKS = -I/opt/local/include /opt/local/lib/libCGAL.* /opt/local/lib/libgmp.*
	OPTFLAGS = -w -stdlib=libc++
endif

#DIR1 = includes

#DIR2 = data_structures

#DIR3 = image_manip

# DIR4 = helpers

# DIR5 = /usr/local/include/eigen3

# INC = $(DIR1) $(DIR2) $(DIR3) $(DIR4) $(DIR5)

# INC_PARAMS = $(foreach d, $(INC), -I$d)

SOURCES = main.cpp mesh.cpp

OBJECT = geommodel.o

all : clean build run

clean :
	rm -rf opengl.o

build :
	$(CXX) $(SOURCES) $(LIBS) $(OPTFLAGS) $(BADLINKS) -o $(OBJECT)

run :
	./geommodel.o

build_dbg :
	$(CXX) $(SOURCES) $(LIBS) $(OPTFLAGS) $(BADLINKS) -o $(OBJECT)  -g

run_dbg :
	$(GDB) geommodel.o

dbg : clean build_dbg run_dbg
