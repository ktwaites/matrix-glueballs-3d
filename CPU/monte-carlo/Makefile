FC = gfortran
FFLAGS ?= -O2
LDLIBS ?= -llapack -lblas

BUILD_DIR = build
COMMON_OBJECTS = $(BUILD_DIR)/matrix-glueball-3d-parameters.o \
                 $(BUILD_DIR)/elementary_functions.o \
                 $(BUILD_DIR)/linalg.o $(BUILD_DIR)/ranlux.o \
                 $(BUILD_DIR)/lattice.o $(BUILD_DIR)/monte_carlo.o

.PHONY: all modules clean

all: save_frozen thermalise convert_csv

modules: $(COMMON_OBJECTS)

save_frozen: $(BUILD_DIR)/frozen_field.o $(COMMON_OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

thermalise: $(BUILD_DIR)/thermalise.o $(COMMON_OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

convert_csv: $(BUILD_DIR)/convert_csv.o $(COMMON_OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/matrix-glueball-3d-parameters.o: matrix-glueball-3d-parameters.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/elementary_functions.o: elementary_functions.f90 $(BUILD_DIR)/matrix-glueball-3d-parameters.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/linalg.o: linalg.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/ranlux.o: ranlux.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/lattice.o: lattice.f90 $(BUILD_DIR)/matrix-glueball-3d-parameters.o $(BUILD_DIR)/elementary_functions.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/monte_carlo.o: monte_carlo.f90 $(BUILD_DIR)/matrix-glueball-3d-parameters.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/ranlux.o $(BUILD_DIR)/lattice.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/frozen_field.o: frozen_field.f90 $(BUILD_DIR)/matrix-glueball-3d-parameters.o $(BUILD_DIR)/lattice.o $(BUILD_DIR)/monte_carlo.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/thermalise.o: thermalise.f90 $(BUILD_DIR)/matrix-glueball-3d-parameters.o $(BUILD_DIR)/lattice.o $(BUILD_DIR)/monte_carlo.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/convert_csv.o: convert_csv.f90 $(BUILD_DIR)/matrix-glueball-3d-parameters.o $(BUILD_DIR)/monte_carlo.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o *.mod save_frozen thermalise convert_csv
