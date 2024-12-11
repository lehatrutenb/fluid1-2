run:
	g++ -DDEBUG="true" -std=c++2b fluid.cpp
	./a.out

runTest:
	g++ -DDEBUG="true" -std=c++2b test.cpp
	./a.out