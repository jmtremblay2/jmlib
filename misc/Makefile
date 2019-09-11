all: install

CPP = g++
CFLAGS = -c -O2 -Wall -Werror -fpic

#C = gcc
#CFLAGS = -c O2 -Wall -Werror -fpic

INSTALLPATH=/usr/local
FP=/usr/local

SRC = $(wildcard *.cpp)
HEAD = $(wildcard *.hpp)
OBJS = $(SRC:.cpp=.o)
 
%.o : %.cpp
	$(CPP) $(CFLAGS) -o $@ $<

mvt_f2c.o: fortran/mvt_f2c.cpp
	$(CPP) -fPIC fortran/mvt_f2c.cpp -c -O3 -o mvt_f2c.o
	
libjm.so: $(OBJS) mvt_f2c.o
	g++ -shared -o libjm.so $(OBJS) mvt_f2c.o \
			/usr/lib/libarmadillo.so /usr/lib/x86_64-linux-gnu/libnlopt.so -lm -lf2c

install: libjm.so
	sudo mkdir -p $(INSTALLPATH)/include/jm
	sudo cp $(HEAD) fortran/mvt_f2c.hpp $(INSTALLPATH)/include/jm/
	sudo cp libjm.so $(INSTALLPATH)/lib/

clean :
	@rm -f *.o *.so *~
