CC=g++
CFLAGS=-c -Wall -std=c++11
LDFLAGS=-lgsl -lgslcblas
SOURCES=mcml_main.cpp mcml_model.cpp mcml_conv.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mcml

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)
