main:
	g++ -O1 main_with_par.cpp -o main -fopenmp
test:
	g++ --std=c++11 vector_operations_tests.cpp -o test -I googletest/googletest/include googletest/lib/libgtest.a -pthread
clean:
	rm -f main test
