all: pagerank.cpp 
	g++ -o pagerank pagerank.cpp

clean: 
	$(RM) pagerank

