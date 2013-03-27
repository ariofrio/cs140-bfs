OBJS=sequential parallel0 parallel1 parallel2

all: $(OBJS) rmat rmat_matlab

parallel0: parallel.cpp
	cilk++ -Dmethod=0  -o $@ $< -lcilkutil

parallel1: parallel.cpp
	cilk++ -Dmethod=1  -o $@ $< -lcilkutil

parallel2: parallel.cpp
	cilk++ -Dmethod=2  -o $@ $< -lcilkutil

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
