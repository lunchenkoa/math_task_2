main: 
	g++ main.cpp iof.cpp functions.cpp -o main
clean:
	rm -fr *.o main ./solution/*.txt
result: main
	./main
	python3 ./result_processing/graphs.py
