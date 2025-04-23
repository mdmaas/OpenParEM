bin_PROGRAMS = OpenParEM2D simplify2D process2D
noinst_PROGRAMS = waveguide

OpenParEM2D_SOURCES = convergence.cpp fem2D.cpp fieldPoints.cpp findPoints.cpp modes.cpp OpenParEM2D.cpp results.cpp eigensolve.c Hsolve.c project.c convergence.hpp fieldPoints.hpp fem2D.hpp fieldPoints.hpp modes.hpp results.hpp eigensolve.h Hsolve.h project.h 

waveguide_SOURCES = waveguide.c waveguide.h

simplify2D_SOURCES = simplify2D.cpp project.c simplify2D.hpp project.h

process2D_SOURCES = process2D.cpp project.c process2D.hpp project.h

check-local:
	../../scripts/OpenParEM2D_test
