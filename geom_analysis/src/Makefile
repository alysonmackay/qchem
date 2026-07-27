CXX	 = g++ 
CXXFLAGS = -std=c++14 -Wall -O2
INCLUDES = -isystem ../eigen 

TARGET = geom 
OBJECTS = geom.o molecule.o 

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp molecule.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean: 
	rm -f $(OBJECTS) $(TARGET) acetaldehyde.out

run: $(TARGET) 
	./$(TARGET) ../input/acetaldehyde.dat >> acetaldehyde.out
