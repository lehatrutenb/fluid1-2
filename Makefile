run:
	g++ -DDEBUG="true" -std=c++2b fluid.cpp
	./a.out

runTest:
	g++ -DDEBUG="true" -std=c++2b test.cpp
	./a.out

runTest2:
	g++ -DDEBUG="true" -std=c++2b -DTYPES="FIXED(10, 10),FIXED(20, 20)" test2.cpp
	./a.out

runTest2E:
	g++ -DDEBUG="true" -std=c++2b -E -DTYPES="FIXED(10, 10),FIXED(20, 20),FLOAT" test2.cpp
	./a.out

runTest3:
	g++ -DDEBUG="true" -std=c++2b test3.cpp
	./a.out