CC = cc
CFLAGS ?= -O2 -g
WARN = -std=c23 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
OPENMP = -fopenmp
OBJ = main.o graph.o steiner.o topology.o realize.o solve.o narrow.o embed.o reader.o

.PHONY: all test clean
all: planar-cut

planar-cut: $(OBJ)
	$(CC) $(LDFLAGS) $(OPENMP) -o $@ $(OBJ) $(LDLIBS)

%.o: %.c cut.h
	$(CC) $(CPPFLAGS) $(WARN) $(CFLAGS) $(OPENMP) -c -o $@ $<

embed-test: tests/embed_test.o embed.o graph.o steiner.o
	$(CC) $(LDFLAGS) $(OPENMP) -o $@ $^ $(LDLIBS)

test: planar-cut embed-test
	./embed-test
	python3 tests/test.py ./planar-cut
	python3 tests/test.py ./planar-cut --klein-marx-only

clean:
	rm -f $(OBJ) tests/embed_test.o planar-cut embed-test
