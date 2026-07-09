#(c) Copyright Ruben Cardenes Almeida 13-9-2002
#
# Sources live in src/, object files are built into build/,
# and executables are produced in the project root.

CC      = gcc
SRC     = src
BUILD   = build
# -std=gnu89 + -Wno-implicit-function-declaration keep this 2004-era C
# building on modern clang/gcc (implicit declarations are otherwise errors).
FLAGS   = -g -O2 -std=gnu89 -Wno-implicit-function-declaration -I$(SRC)

# png_write.c instantiates the stb_image single-header libraries, which require
# C99 or later, so it is compiled with its own flags.
PNGFLAGS = -g -O2 -std=gnu11 -I$(SRC)

# Executables listed in README.md
EXECUTABLES = \
	thickness2D \
	thickness2D_knee \
	thickness3D \
	laplace2D \
	laplace3D \
	poisson2D \
	compute_boundary2D \
	compute_boundary_knee3D \
	compute_boundary_cortex3D \
	generar_phantom3d_elipsoid

all: $(EXECUTABLES)

## Object files (built into build/)
$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/thickness2D.o: $(SRC)/thickness2D.c $(SRC)/thickness2D.h | $(BUILD)
	$(CC) $(FLAGS) -c $(SRC)/thickness2D.c -o $@

$(BUILD)/thickness3D.o: $(SRC)/thickness3D.c $(SRC)/thickness3D.h | $(BUILD)
	$(CC) $(FLAGS) -c $(SRC)/thickness3D.c -o $@

$(BUILD)/laplace2D.o: $(SRC)/laplace2D.c $(SRC)/laplace2D.h | $(BUILD)
	$(CC) $(FLAGS) -c $(SRC)/laplace2D.c -o $@

$(BUILD)/laplace3D.o: $(SRC)/laplace3D.c $(SRC)/laplace3D.h | $(BUILD)
	$(CC) $(FLAGS) -c $(SRC)/laplace3D.c -o $@

$(BUILD)/poisson2D.o: $(SRC)/poisson2D.c $(SRC)/poisson2D.h | $(BUILD)
	$(CC) $(FLAGS) -c $(SRC)/poisson2D.c -o $@

$(BUILD)/io.o: $(SRC)/io.c $(SRC)/io.h | $(BUILD)
	$(CC) $(FLAGS) -c $(SRC)/io.c -o $@

$(BUILD)/png_write.o: $(SRC)/png_write.c $(SRC)/png_write.h \
                      $(SRC)/stb_image/stb_image.h $(SRC)/stb_image/stb_image_write.h | $(BUILD)
	$(CC) $(PNGFLAGS) -c $(SRC)/png_write.c -o $@

$(BUILD)/DToptimo3d.o: $(SRC)/DToptimo3d.c $(SRC)/DToptimo3d.h | $(BUILD)
	$(CC) $(FLAGS) -c $(SRC)/DToptimo3d.c -o $@

## Executables (produced in the project root)

thickness2D: $(SRC)/thickness2Dmain.c $(BUILD)/thickness2D.o $(BUILD)/laplace2D.o $(BUILD)/png_write.o
	$(CC) $(FLAGS) $(SRC)/thickness2Dmain.c $(BUILD)/laplace2D.o $(BUILD)/thickness2D.o $(BUILD)/png_write.o -lm -o $@

thickness2D_knee: $(SRC)/thickness2Dmain_knee.c $(BUILD)/thickness2D.o $(BUILD)/laplace2D.o $(BUILD)/png_write.o
	$(CC) $(FLAGS) $(SRC)/thickness2Dmain_knee.c $(BUILD)/laplace2D.o $(BUILD)/thickness2D.o $(BUILD)/png_write.o -lm -o $@

thickness3D: $(SRC)/thickness3Dmain.c $(BUILD)/laplace3D.o $(BUILD)/thickness3D.o $(BUILD)/io.o $(BUILD)/DToptimo3d.o
	$(CC) $(FLAGS) $(SRC)/thickness3Dmain.c $(BUILD)/laplace3D.o $(BUILD)/thickness3D.o $(BUILD)/io.o $(BUILD)/DToptimo3d.o -lm -o $@

laplace2D: $(SRC)/laplace2Dmain.c $(BUILD)/laplace2D.o $(BUILD)/png_write.o
	$(CC) $(FLAGS) $(SRC)/laplace2Dmain.c $(BUILD)/laplace2D.o $(BUILD)/png_write.o -lm -o $@

laplace3D: $(SRC)/laplace3Dmain.c $(BUILD)/laplace3D.o
	$(CC) $(FLAGS) $(SRC)/laplace3Dmain.c $(BUILD)/laplace3D.o -lm -o $@

poisson2D: $(SRC)/poisson2Dmain.c $(BUILD)/poisson2D.o $(BUILD)/laplace2D.o
	$(CC) $(FLAGS) $(SRC)/poisson2Dmain.c $(BUILD)/poisson2D.o $(BUILD)/laplace2D.o -lm -o $@

compute_boundary2D: $(SRC)/compute_boundary2D.c $(BUILD)/io.o $(BUILD)/laplace2D.o
	$(CC) $(FLAGS) $(SRC)/compute_boundary2D.c $(BUILD)/laplace2D.o $(BUILD)/io.o -lm -o $@

compute_boundary_knee3D: $(SRC)/compute_boundary_knee3D.c $(BUILD)/io.o $(BUILD)/laplace3D.o
	$(CC) $(FLAGS) $(SRC)/compute_boundary_knee3D.c $(BUILD)/laplace3D.o $(BUILD)/io.o -lm -o $@

compute_boundary_cortex3D: $(SRC)/compute_boundary_cortex3D.c $(BUILD)/io.o $(BUILD)/laplace3D.o
	$(CC) $(FLAGS) $(SRC)/compute_boundary_cortex3D.c $(BUILD)/laplace3D.o $(BUILD)/io.o -lm -o $@

generar_phantom3d_elipsoid: $(SRC)/generar_phantom3d_elipsoid.c
	$(CC) $(FLAGS) $(SRC)/generar_phantom3d_elipsoid.c -o $@

clean:
	rm -f $(EXECUTABLES)
	rm -rf $(BUILD)
	rm -rf *.dSYM
	rm -f *~
	rm -f *.volf
	rm -f *.vols
	rm -f *.flt
	rm -f core
