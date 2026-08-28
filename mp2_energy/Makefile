CXX      = g++
CXXFLAGS = -std=c++14 -Wall -O2
# c++14 req for old eigen version
INCLUDES = -isystem ../eigen/ -I../hartreefock_scf/ -I../tools/

TARGET  = main
OBJECTS = main.o mp2.o ../hartreefock_scf/hf_scf.o ../tools/orb_transform.o

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean:
	rm -f $(OBJECTS) $(TARGET) mp2_h2o.out

test: $(TARGET) 
	./$(TARGET) h2o STO-3G >> mp2_h2o.out 
