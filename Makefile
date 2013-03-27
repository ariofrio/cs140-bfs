CXX=cilk++
CXXFLAGS=-lcilkutil
OBJS=sequential parallel0 parallel1 parallel2 parallel3

all: $(OBJS) rmat rmat_matlab

parallel0: parallel.cpp
	$(CXX) -Dmethod=0 -o $@ $< $(CXXFLAGS)

parallel1: parallel.cpp
	$(CXX) -Dmethod=1 -o $@ $< $(CXXFLAGS)

parallel2: parallel.cpp
	$(CXX) -Dmethod=2 -o $@ $< $(CXXFLAGS)

parallel3: parallel.cpp
	$(CXX) -Dmethod=3 -o $@ $< $(CXXFLAGS)

rmat: rmat.c
	cc -o rmat rmat.c

rmat_matlab: rmat.c
	cc -DMATLAB -o rmat_matlab rmat.c

clean:
	rm -f $(OBJS) rmat rmat_matlab

.DEFAULT_GOAL=test
test: all
	@./test.sh

sync:
	while true; do rsync --archive . triton:bfs; sleep 0.1; done
