cgtest: main.o CGSolve.o InnerProduct.o
	g++ main.o CGSolve.o InnerProduct.o -lblas -llapack

main.o: main.cpp
	g++ -c main.cpp

CGSolve.o: CGSolve.cpp
	g++ -c CGSolve.cpp

InnerProduct.o: InnerProduct.cpp
	g++ -c InnerProduct.cpp
