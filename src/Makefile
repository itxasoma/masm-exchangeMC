FC      = gfortran
FFLAGS  = -O2 -g -Wall -Wextra -fcheck=all -fbacktrace
CPPFLAGS = -I./r1279
TARGET  = exchange_mc          # Parts 1, 2, 4, 5. time-series output
TARGET3 = exchange_mc_part3    # Part 3. histogram output
TARGET6 = exchange_mc_part6    # Part 6. time-series output for NT=1
PYTHON  = python3

SHARED_OBJ = parameters.o lattice.o rng_wrapper.o bonds.o exchangeMC.o r1279.o ran2.o

all: $(TARGET) $(TARGET3) $(TARGET6)

$(TARGET): $(SHARED_OBJ) main.o
	$(FC) $(FFLAGS) -o $@ $^

$(TARGET3): $(SHARED_OBJ) main3.o
	$(FC) $(FFLAGS) -o $@ $^

$(TARGET6): $(SHARED_OBJ) main6.o
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

main3.o: main3.f90 parameters.o lattice.o bonds.o exchangeMC.o rng_wrapper.o
	$(FC) $(FFLAGS) -c $<

main6.o: main6.f90 parameters.o lattice.o bonds.o exchangeMC.o rng_wrapper.o
	$(FC) $(FFLAGS) -c $<

ran2.o: r1279/ran2.f
	$(FC) $(FFLAGS) -c $< -o $@

r1279.o: r1279/r1279.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@

# To run the code:
run: $(TARGET)
	./$(TARGET) ../inputs/part1.in
	mv *.dat ../results/

run2: $(TARGET)
	./$(TARGET) ../inputs/part2.in
	mv *.dat ../results/

run3: $(TARGET3)
	./$(TARGET3) ../inputs/part3.in
	mv *.dat ../results/

run6: $(TARGET6)
	mkdir -p ../results/part6
	./$(TARGET6) ../inputs/part6.in
	mv timeseries_part6*.dat ../results/part6/

# Figures
figs2:
	mkdir -p ../results ../figures
	$(PYTHON) plots2.py

figs3:
	mkdir -p ../results ../figures
	$(PYTHON) plots3.py

figs4:
	mkdir -p ../figures
	$(PYTHON) plots4.py

figs5:
	mkdir -p ../figures
	$(PYTHON) plots5.py

figs6:
	mkdir -p ../figures
	$(PYTHON) plots6.py

figs7:
	mkdir -p ../figures
	$(PYTHON) plots7.py

# Clean all
clean:
	rm -f *.o *.mod $(TARGET) $(TARGET3) $(TARGET6) *.x fort.*


# Binnings
BIN_TARGET = main_binning.x

binning.o: binning.f90
	$(FC) $(FFLAGS) -c $<

main_binning.o: main_binning.f90 binning.o
	$(FC) $(FFLAGS) -c $<

$(BIN_TARGET): binning.o main_binning.o
	$(FC) $(FFLAGS) -o $@ $^

bin2: $(BIN_TARGET)
	./$(BIN_TARGET)
	

BIN4_TARGET = main_binning_part4.x

main_binning_part4.o: main_binning_part4.f90 binning.o
	$(FC) $(FFLAGS) -c $<

$(BIN4_TARGET): binning.o main_binning_part4.o
	$(FC) $(FFLAGS) -o $@ $^

bin4: $(BIN4_TARGET)
	@for f in ../results/part4/timeseries_part4_s*.dat; do \
	  tag=$$(basename $$f .dat | sed 's/timeseries_//'); \
	  ./$(BIN4_TARGET) $$f ../results/part4/summary_$${tag}.dat; \
	done


# Autocorrelation
AUTOCORR_TARGET = main_autocorr_part7.x

main_autocorr_part7.o: main_autocorr_part7.f90
	$(FC) $(FFLAGS) -c $<

$(AUTOCORR_TARGET): main_autocorr_part7.o
	$(FC) $(FFLAGS) -o $@ $^

autocorr7: $(AUTOCORR_TARGET)