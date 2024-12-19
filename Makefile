run:
	g++ -DDEBUG="true" -DTYPES="FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT" -DSIZES="S(1920,1080),S(36,84)" -std=c++2b fluid.cpp
	./a.out --p-type="FLOAT" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FIXED(20, 10)" --n=36 --m=84 --field-file="field"

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