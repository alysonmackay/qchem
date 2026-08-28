CXX      = g++
CXXFLAGS = -std=c++14 -Wall -O2
INCLUDES = -isystem ../eigen/ -I../hartreefock_scf/ -I../mp2_energy/

TARGET  = main
OBJECTS = main.o ccsd.o ccsd_setup.o cc_int.o ../mp2_energy/mp2.o ../hartreefock_scf/hf_scf.o

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run
clean:
	rm -f $(OBJECTS) $(TARGET) ccsd_h2o.out

run: $(TARGET) 
	./$(TARGET) h2o STO-3G >> ccsd_h2o.out 
