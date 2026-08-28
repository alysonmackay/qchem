CXX      = g++
CXXFLAGS = -std=c++14 -Wall -O2
INCLUDES = -isystem ../../eigen/

TARGET  = vib_analysis
OBJECTS = vib_analysis.o molecule.o

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp molecule.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean:
	rm -f $(OBJECTS) $(TARGET) h2o_vib_freq.out

run: $(TARGET) 
	./$(TARGET) ../input/h2o_geom.txt 0 ../input/h2o_hessian.txt >> h2o_vib_freq.out 
