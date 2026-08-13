CXX      = g++
CXXFLAGS = -std=c++14 -Wall -O2
INCLUDES = -isystem ../eigen/ -I../hartreefock_scf/

TARGET  = main
OBJECTS = main.o mp2.o ../hartreefock_scf/hf_scf.o

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean:
	rm -f $(OBJECTS) $(TARGET) mp2_h2o.out

run: $(TARGET) 
	./$(TARGET) h2o STO-3G >> mp2_h2o.out 
