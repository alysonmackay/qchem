CXX      = g++
CXXFLAGS = -std=c++14 -Wall -O2
INCLUDES = -isystem ../eigen/ -I../hartreefock_scf/

TARGET  = mp2_energy
OBJECTS = mp2_energy.o ../hartreefock_scf/hf.o

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean:
	rm -f $(OBJECTS) $(TARGET) mp2.out

run: $(TARGET) 
	./$(TARGET) h2o STO-3G >> mp2.out 
