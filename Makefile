# Source files
SOURCES = mod_global_variables.f90 mod_global_functions.f90 mod_calculate_rates.f90 dvode_f90_m.f90 \
          mod_read_model.f90 mod_read_rate06.f90 mod_save_results.f90 mod_run_dvode.f90 chem_rate06_dvode.f90

# Object files (derived from source files)
OBJECTS = $(SOURCES:.f90=.o)

# Compiler and flags
#FC = ifort
 FC = gfortran

#FLAGS = -g -pg -O2 -fprotect-parens -fp-model=strict -march=native -module . -diag-disable=10448  # ifort options

# Chemistry/physics modules: O3 + native tuning + real*8 promotion
FLAGS = -O3 -march=native -funroll-loops -fdefault-real-8 -ffree-line-length-512 -J .

# DVODE gets conservative flags: no -O3, no -ffast-math, strict FP behavior
DVODE_FLAGS = -O2 -march=native -fdefault-real-8 -ffree-line-length-512 -J .

# Debug build: swap FLAGS -> DEBUG_FLAGS and recompile
DEBUG_FLAGS = -g -pg -O0 -fbacktrace -fcheck=all -ffpe-trap=invalid,overflow,zero \
              -fdefault-real-8 -ffree-line-length-512 -J . -Wall -Wextra

# Target executable
TARGET = monaco

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) -o $@ $(OBJECTS)      # TAB required here
	@echo "Make complete"                # TAB required here

# Compilation rules with dependencies

mod_global_variables.o: mod_global_variables.f90
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_global_functions.o: mod_global_functions.f90 mod_global_variables.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_calculate_rates.o: mod_calculate_rates.f90 mod_global_variables.o mod_global_functions.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

# DVODE compiled separately with conservative flags
dvode_f90_m.o: dvode_f90_m.f90
	$(FC) $(DVODE_FLAGS) -c $< -o $@     # TAB required here

mod_read_rate06.o: mod_read_rate06.f90 mod_global_variables.o mod_global_functions.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_read_model.o: mod_read_model.f90 mod_global_variables.o mod_global_functions.o mod_read_rate06.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_save_results.o: mod_save_results.f90 mod_global_variables.o mod_global_functions.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_run_dvode.o: mod_run_dvode.f90 mod_global_variables.o mod_global_functions.o mod_calculate_rates.o dvode_f90_m.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

chem_rate06_dvode.o: chem_rate06_dvode.f90 mod_global_variables.o mod_global_functions.o mod_calculate_rates.o \
                     dvode_f90_m.o mod_read_model.o mod_read_rate06.o mod_save_results.o mod_run_dvode.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

# Clean up
clean:
	rm -rf *~ *.o *.mod ab csv fort.* results* analytics*  # TAB required here

# Phony targets
.PHONY: all clean
