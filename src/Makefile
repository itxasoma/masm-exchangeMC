FC = gfortran
FFLAGS = -O2 -g -Wall -Wextra -fcheck=all -fbacktrace
CPPFLAGS = -I./r1279
TARGET = exchange_mc
PYTHON = python3

OBJ = parameters.o lattice.o rng_wrapper.o bonds.o exchangeMC.o main.o r1279.o ran2.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

parameters.o: parameters.f90
	$(FC) $(FFLAGS) -c $<

lattice.o: lattice.f90 parameters.o
	$(FC) $(FFLAGS) -c $<

rng_wrapper.o: rng_wrapper.f90
	$(FC) $(FFLAGS) -c $<

bonds.o: bonds.f90 parameters.o lattice.o rng_wrapper.o
	$(FC) $(FFLAGS) -c $<

exchangeMC.o: exchangeMC.f90 parameters.o lattice.o bonds.o rng_wrapper.o
	$(FC) $(FFLAGS) -c $<

main.o: main.f90 parameters.o lattice.o bonds.o exchangeMC.o rng_wrapper.o
	$(FC) $(FFLAGS) -c $<

ran2.o: r1279/ran2.f
	$(FC) $(FFLAGS) -c $< -o $@

r1279.o: r1279/r1279.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET) ../inputs/part1.in
	mv *.dat ../results/

run2: $(TARGET)
	./$(TARGET) ../inputs/part2.in
	mv *.dat ../results/

figures:
	mkdir -p ../results ../figures
	python3 plots.py

clean:
	rm -f *.o *.mod $(TARGET) *.x fort.*


BIN_TARGET = main_binning.x

binning.o: binning.f90
	$(FC) $(FFLAGS) -c $<

main_binning.o: main_binning.f90 binning.o
	$(FC) $(FFLAGS) -c $<

$(BIN_TARGET): binning.o main_binning.o
	$(FC) $(FFLAGS) -o $@ $^

binning: $(BIN_TARGET)
	./$(BIN_TARGET)