CXX      = g++
CXXFLAGS = -std=c++14 -Wall -O2
INCLUDES = -isystem ../eigen/ -I../tools/ -I../hartreefock_scf/

TARGET  = main
OBJECTS = main.o triples.o ../tools/orb_transform.o ccsd.o ccsd_setup.o cc_int.o ../hartreefock_scf/hf_scf.o 

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean:
	rm -f $(OBJECTS) $(TARGET) ccsd_h2o.out

test: $(TARGET) 
	./$(TARGET) h2o STO-3G >> ccsd_h2o.out 
