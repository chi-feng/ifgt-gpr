bin/ifgt : src/ifgt.cpp src/IFGT.hh src/linalg.hh src/rng.hh
	g++ src/ifgt.cpp -o bin/ifgt -Wall -Werror -pedantic -O3