g++ -O2 -pthread -std=c++11 -I./include -I./wsserver -I/usr/include/eigen3/ -L./bin -c dynamics.cpp 
ar rcs ./bin/libdynamics.a dynamics.o
rm *.o