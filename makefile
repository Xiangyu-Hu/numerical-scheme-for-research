OBJECTS = fluid.o init.o main.o output.o solver.o glbfunc.o

COMPILER = g++

OPT = -O2  -Wall


program: $(OBJECTS)
	$(COMPILER) $(OPT) $(OBJECTS) -o run_PGNSV
fluid.o: fluid.cpp glbcls.h glbfunc.h
	$(COMPILER) $(OPT) -c $(INCLUDEPATH) fluid.cpp 
init.o: init.cpp glbcls.h glbfunc.h
	$(COMPILER) $(OPT) -c init.cpp 
main.o: main.cpp glbcls.h glbfunc.h
	$(COMPILER) $(OPT) -c main.cpp 
utput.o: output.cpp glbcls.h glbfunc.h
	$(COMPILER) $(OPT) -c $(INCLUDEPATH) output.cpp 
solver.o: solver.cpp glbcls.h glbfunc.h
	$(COMPILER) $(OPT) -c solver.cpp 
glbfunc.o: glbfunc.cpp glbfunc.h
	$(COMPILER) $(OPT) -c glbfunc.cpp 

.PHONY: clean
clean:
	@rm -f run_PGNSV $(OBJECTS)


# ============================================================
#        SOME EXPLANATIONS TO THIS MAKEFILE
#        ----------------------------------
# --> Define variables with Var=... and link to these with $()
# --> Define dependancies via "target: dependant objects".
# 	Watch for: Execution lines below with tabulator!
# --> Target is updated only if dependancies are changed;
# 	therefore "clean", which is independant, must be 	
#	forced to always be executed (with "make clean") 
# 	via the ".PHONY: clean" command.
# --> Explanation of options in compiler command:
# 	-Wall	outputs warnings
# 	-O3	activates (rigid) optimization
# 	-g0	activates debugging (not recommended with -O3)
#	-openmp	activates parallelization
#	-o 	creates output-file
# 	-c	creates .o-file to be used later
#	@	suppresses output on display
#	-f	ignores error messages
# ============================================================
