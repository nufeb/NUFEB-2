SRCS := $(wildcard *.cpp)
OBJS := $(addsuffix .o, $(basename $(SRCS)))
GUI_EXE = minimal

all: $(GUI_EXE)

#GeneNet creation Tool Chain
$(GUI_EXE): $(SRCS)
	 g++ minimal.cpp `wx-config-3.0 --cxxflags --libs` -o minimal

.PHONY: clean

clean:
	rm -rf $(OBJS) $(GUI_EXE) *~
