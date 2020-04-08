CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/gsl/
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lm
SF_OBJECTS	= simulation_run_eff.o utilities_comp.o io.o

eff:	$(SF_OBJECTS)
	$(CC) $(CC_FLAGS) $(SF_OBJECTS) -o comp_sim_eff $(LD_FLAGS)
utilities_comp.o: utilities_comp.cpp
	$(CC) $(CC_FLAGS) -c utilities_comp.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
simulation_run_eff.o: simulation_run_eff.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_eff.cpp

clean:
	\rm *.o *~ 

