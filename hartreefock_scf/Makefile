CXX      = g++
CXXFLAGS = -std=c++23 -Wall -O2
#INCLUDES = -isystem ../../eigen/

TARGET  = hf_scf
OBJECTS = hf_scf.o 

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean:
	rm -f $(OBJECTS) $(TARGET) hf_scf.out

run: $(TARGET) 
	./$(TARGET) input/h2o/STO-3G/enuc.dat >> hf_scf.out 
