main: 
	 g++ de_allocate.cpp eq_sol_id_gas.cpp main.cpp -o main
clean:
	rm -fr *.o main
result: main
	./main
	python3 ./result_processing/graphs.py
