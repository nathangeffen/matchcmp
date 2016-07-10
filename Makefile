CXX = g++

CXXFLAGS = -Wall -std=c++11
DEVFLAGS  = -g -rdynamic
RELFLAGS = -O3
LDFLAGS =

# the build target executable:
SOURCES = partners.cc
OBJECTS = $(SOURCES:.cc=.o)
DEPEND =  $(OBJECTS:%.o=.%.d)
EXECUTABLE = partners

all: $(SOURCES) $(EXECUTABLE)-dev

$(EXECUTABLE)-dev: $(OBJECTS)
	$(CXX) $(DEVFLAGS) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)-dev

%.o: %.cc
	$(CXX) -c $(DEVFLAGS) $(CXXFLAGS)  -MD -MP -MF .${@:.o=.d} $< -o $@

release: clean
	$(CXX) $(RELFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $(EXECUTABLE)-rel $(SOURCES)

clean:
	rm -f $(EXECUTABLE)-rel $(EXECUTABLE)-dev *.o

-include $(DEPEND)
