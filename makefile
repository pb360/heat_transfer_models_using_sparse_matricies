main: main.cpp  CGSolver.cpp CGSolver.hpp COO2CSR.cpp COO2CSR.hpp heat.cpp heat.hpp matvecops.cpp matvecops.hpp sparse.cpp sparse.hpp 
	g++ -std=c++11 -Wall -Wconversion -Wextra -Wpedantic -Wunused-parameter -Wunused-parameter -Wuninitialized -o main main.cpp sparse.cpp COO2CSR.cpp CGSolver.cpp matvecops.cpp

.PHONY: clean
clean:
	$(RM) main *.o