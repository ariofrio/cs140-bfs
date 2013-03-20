OBJS=sequential parallel

all: $(OBJS)

parallel: parallel.cpp
	cilk++ -fcilkscreen -o $@ $<

clean:
	rm -f $(OBJS)

.DEFAULT_GOAL=test
test: all
	./test.sh

sync:
	while true; do rsync --archive . triton:bfs; sleep 0.1; done
