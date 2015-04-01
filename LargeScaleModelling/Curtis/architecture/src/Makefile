SRCS := $(wildcard *.cpp)
OBJS := $(addsuffix .o, $(basename $(SRCS)))
NUFEB_EXE = ../bin/nufeb

all: $(NUFEB_EXE)

#GeneNet creation Tool Chain
$(NUFEB_EXE): $(SRCS) $(OBJS)
	g++ -o nufeb $(OBJS)
	mv nufeb ../bin/

%.o: %.cpp
	g++ -O3 -Wall -c -fmessage-length=0 -o$@ $<

.PHONY: clean

clean:
	rm -rf $(OBJS) $(NUFEB_EXE) *~
