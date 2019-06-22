##########################################################
############### MAKEFILE FOR THE MINFINDER ###############
OBJS = base.o problem.o main.o getoptions.o doublepop.o cobyla.o tolmin.o

include ../Makefile.inc

all: $(OBJS)  make_program

%.o: %.cc
	$(CXX)  $(CFLAGS) -I../include -c $< $(INC) -o $@

make_program: make_program.cc
	$(CXX) make_program.cc -DCC=\"$(CC)\" -DCXX=\"$(CXX)\" -DF77=\"$(F77)\" -DF77FLAGS=\"$(F77FLAGS)\" -DROOTDIR=\"$(ROOTDIR)\" -DLINK=\"$(LINK)\" -DINCL=\"$(INCL)\" -o make_program
	mv make_program ../bin

clean:
	rm -f *~ *.o   make_program
##########################################################
##########################################################
