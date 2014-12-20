CC=g++
CFLAGS=-c -Wall -std=c++11
LDFLAGS=
SOURCES=mcml_main.cpp mcml_model.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mcml

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)
