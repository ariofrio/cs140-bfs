OBJS=sequential parallel0 parallel1

all: $(OBJS)

parallel0: parallel.cpp
	cilk++ -Dmethod=0 -fcilkscreen -o $@ $< -lcilkutil

parallel1: parallel.cpp
	cilk++ -Dmethod=1 -fcilkscreen -o $@ $< -lcilkutil

parallel2: parallel.cpp
	cilk++ -Dmethod=2 -fcilkscreen -o $@ $< -lcilkutil

clean:
	rm -f $(OBJS)

.DEFAULT_GOAL=test
test: all
	./test.sh

sync:
	while true; do rsync --archive . triton:bfs; sleep 0.1; done
