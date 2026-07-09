# (c) Copyright, Rub�n C�rdenes Almeida. March, 2004
# Thikness algorithms (2D and 3D) using a Laplacian scalar field 
# (solving PDE's proposed by Yezzi) and using ordered propagation

## Build
```sh
make
```
Builds all executables into the project root (objects go under `build/`).

## Image formats (2D tools)
The 2D tools (`thickness2D`, `thickness2D_knee`, `laplace2D`) read their input
domain from an **8-bit grayscale PNG** and take the image dimensions from the
file, so `nrows`/`ncols` are no longer passed on the command line. Output is a
PNG when the output filename ends in `.png` (grayscale, or an RGB color map with
`-c 1`/`-c 2`), or a raw float file for any other extension. The `.chr` domains
in `data/` have equivalent `.png` versions (see `data/*.png`). The 3D tools
still use raw `.vols`/`.vol` volumes.

-----------------------
# thickness2D

Computes the 2D thickness of a band-shaped region (for example the cerebral
cortex) between its inner and outer boundaries. It solves the Laplace equation
over the band to build a smooth harmonic field that goes from 0 on one boundary
to 1 on the other, and then integrates the normalized gradient of that field
(the Yezzi–Prince PDE method) to measure, at every pixel, the distance from one
boundary to the other. The two boundaries are derived automatically from a
tissue segmentation: the interior is given by the white-matter label (`--lw`)
and the band itself by the cortex label (`--lc`).

## Usage 
```
thickness2D [options] input2D.png output2D.png
              -d (debug mode)
              -n iterations_thickness (10)
              -i iterations_laplace (100)
              -w (swap bytes of input) 
              -r (reverse) 
              -s (sum) 
              -l lambda (0.5)
              -m compute and show mean thickness
              -c color output (0: gray, 1: red-blue, 2: random)
              --hx hy (1) 
              --hy hx (1) 
              --lw white matter label (3)
              --lc cortex label (2)
              --DT (compute thickness with Euclidean DT)
              --streamlines (show only the streamlines)


Input: 8-bit grayscale PNG defining the domain, labeled as:
       "lc" in the band
       "lw" in the interior
        whatever different label, in the outer
       The image dimensions (rows, cols) are taken from the PNG, so nrows/ncols
       are no longer passed on the command line.
Output: if the output filename ends in ".png", a PNG image of the thickness
        map (grayscale, or an RGB color map when -c is 1 or 2; thickness
        normalized to 0-255, background = black). Any other extension writes a
        raw float image (no header), as before.
```

## options: 
-d debug option: gives more text information, and saves images for gradients 
   (x,y), laplacian field, and domain
-w swap bytes after reading input file
-r reverse: compute the thickness from 1 to 0 (without -r thickness is 
   computed from 0 to 1)
-s compute thickness in both directions and sum both results and resting 
   one unit
-k knee mode: computation of boundaries for the knee (experimental)
-n iterations thickness: number of iterations for thickness 
   computation (default 10)
-i number of iterations for laplace computation (default 100) 
-l lamba convergence parameter used in laplace computation 
   (0.5 default value is ok) 
-m compute and show mean and standard deviation of thickness 
--lw label of the interior (white matter, default 3) 
--lc label of the band (cortex, default 2)
-c color output for PNG: 0 grayscale (default), 1 red-blue (thin=red, thick=blue),
   2 random colors. Only affects PNG output.
--hx --hy pixel size parameters (default: 1,1)
--DT compute thickness using Euclidean DT (then laplacian is not used)
--streamlines show only the streamlines in the output file (only for testing)

## Examples 2D:
synthetic non uniform ring
```sh
./thickness2D -n 20 -i 100 --hx 1 --hy 1 --lw 3 --lc 2 data/domain_anillo_256.png output_thickness_2d.png
```
synthetic non uniform ring with a high saliency, colored red-blue
```sh
./thickness2D -n 20 -i 100 -c 1 --hx 1 --hy 1 --lw 3 --lc 2 data/domain_anillo_256_2.png output_thickness_2d.png
```
donut: circle inside a circle
```sh
./thickness2D  -n 10 -i 100 --hx 1 --hy 1 --lw 3 --lc 2 data/domain_donut.png output_thickness_2d_donut.png
```
circle inside an ellipse 
```sh
./thickness2D -n 10 -i 100 --hx 1 --hy 1 --lw 3 --lc 2 data/domain_elipse.png output_thickness_2d_elipse.png
```
square inside a square
```sh
./thickness2D -n 10 -i 200 --hx 1 --hy 1 --lw 3 --lc 2 data/domain_cuadrado.png output_thickness_cuadrado.png
```
baby brain
```sh
./thickness2D -n 10 -i 200 --hx 1 --hy 1 --lw 3 --lc 2 $DATA_DIR/baby_brain/segwmgm/A58.png output_thickness.png
```
------------------------------
# thickness2D_knee

Uses the same Laplacian thickness algorithm as `thickness2D`, but for inputs
whose boundaries are already explicitly labeled instead of being derived from a
tissue segmentation. `thickness2D` builds the inner/outer boundaries of the
band itself from the white-matter and cortex labels (via
`compute_boundary_cortex2D`); `thickness2D_knee` skips that step and expects the
domain to already encode the two boundaries directly (as produced by
`compute_boundary2D`), which suits knee-cartilage data where there is no
white-matter interior from which to infer the boundaries. It uses a mostly
positional command line and a fixed band label of 2.

## Usage
```
thickness2D_knee [options] domain2D.png output2D.png iterations_laplace iterations hx hy [reverse] [lambda]
              -c color output (0: gray, 1: red-blue, 2: random)
```
Input: 8-bit grayscale PNG; the image dimensions are read from the file (so
       max1/max2 are no longer passed on the command line).
Output: filename ending in ".png" writes a PNG of the thickness map (grayscale,
        or RGB when -c is 1 or 2; background = black); any other extension writes
        a raw float image, as before.

## Examples
```sh
./thickness2D_knee $DATA_DIR/knee-images/thickness/010/input.022 $DATA_DIR/knee-images/thickness/010/output2D.png 200 10 0.234375 0.234375 1
./thickness2D_knee -c 1 data/domain_knee_512.png output_thickness_knee.png 200 10 1 1 0 
./thickness2D_knee /home/ruben/data/knee-images/thickness/010/input.011 /home/ruben/data/knee-images/thickness/010/output2D.png 200 10 1 1 0 
```

------------------------------
# thickness3D 

Volumetric (3D) version of `thickness2D`. Given a segmented volume, it computes
the inner and outer boundaries of the cortex band from the white-matter (`--lw`)
and cortex (`--lc`) labels, solves the 3D Laplace equation between them, and
integrates the normalized gradient to obtain the thickness of the band at every
voxel. Input is either a single `.vols` unsigned-short volume or a numbered
slice series; output is a float volume of per-voxel boundary-to-boundary
distances.

Usage: 
```sh
thickness3D  [options] segmented.vols/segmented_prefix output3D.volf nrows ncols nslices
              -d (debug mode)
              -n iterations_thickness (10)
              -i iterations_laplace (100)
              -w (swapbytes)
              -r (reverse)
              -s (sum)
              -k (knee mode)
              -l lambda (0.5)
              -m compute and show mean thickness
              --lw white matter label
              --lc cortex label
              --hx hy (1)
              --hy hx (1)
              --hz hz (1)
              --DT (compute thickness with Euclidean DT)

Output: volume float file, representing the distances from one boundary to the other according to a metric defined by the laplacian field
Input: unsigned short volume file (with extension vols), or set of files numbered from 001 to nslices 
       segmented data for the cortex white matter and the rest of classes. 
       The input should be labeled as:
          "lc" in the band to compute the thickness (cortex)
          "lw" in the inner (white matter)
           whatever different label, in the outer (outside the brain)

Output: volume float file, representing the distances from one boundary to the other according to a metric defined by the laplacian field
```

## options: 
-d debug option: gives more text information, and saves volumens for gradients 
   (x,y,z), laplacian field, and domain
-w swap bytes after reading input files
-r reverse: compute the thickness from 1 to 0 (without -r thickness is 
   computed from 0 to 1)
-s compute thickness in both directions and sum both results and resting 
   one unit
-k knee mode: computation of boundaries for the knee (experimental)
-n iterations thickness: number of iterations for thickness computation 
   (default 10)
-i number of iterations for laplace computation (default 100) 
-l lamba convergence parameter used in laplace computation (0.5 default 
   value is ok) 
-m compute and show mean and standard deviation of thickness 
--lw label of the interior (white matter, default 3) 
--lc label of the band (cortex, default 2)
--hx --hy --hz voxel size parameters (default: 1,1,1)
--DT compute thickness using Euclidean DT (then laplacian is not used)

## Examples
```sh
DATA_DIR=/home/common/data/raw/ms/warfieldr21dtmri/1.2.840.113619.2.136.1762888421.2017.1078752767.808.UID
./thickness3D_cortex -m -w -n 10 -i 100 --lw 3 --lc 2 --hx 0.9375 --hy 0.9375 --hz 1.3 $DATA_DIR/brain_124/Output $DATA_DIR/brain_124/Output_thickness.volf 256 256 124 

DATA_DIR=/projects/hpc/active/fernandes/als-neonate-segmented/07002_FG/30th_week/zillesgi/segwmgm
./thickness3D_cortex -m -w -n 10 -i 100 --lw 3 --lc 2 --hx 0.703125 --hy 0.703125 --hz 1.5 $DATA_DIR/baby_brain/segwmgm/I $DATA_DIR/baby_brain/segwmgm/Thickness_zill.volf 256 256 124 

DATA_DIR=/projects/hpc/active/reith/pvl/NMR032_scott_johanna/measure/cerebrum_and_white_matter
./thickness3D_cortex -m -w -n 10 -i 200 --lw 3 --lc 2 --hx 0.703125 --hy 0.703125 --hz 1.5 $DATA_DIR/baby_brain/pvl/cerebrum_and_white_matter $DATA_DIR/baby_brain/pvl/Thickness.volf 256 256 110
```

Datos brainweb:
```sh
./thickness3D -m -n 20 -i 500 --lw 3 --lc 2 ~/datos/datos_meri/knn/out_rf0n3.vols out.volf 187 161 161
```
Caja
```sh
./thickness3D -m -n 20 -i 200 --lw 3 --lc 2 data/input_caja3d.vols Thickness.volf 80 80 80
```
Elipsoid
```sh
./thickness3D -m -n 10 -i 100 --lw 3 --lc 2 data/phantom_elipsoid.vols Thickness.volf 80 80 80
```
Sphere
```sh
./thickness3D -m -n 10 -i 100 --lw 3 --lc 2 data/phantom_sphere.vols Thickness.volf 80 80 80
```
--------------------------------
# thickness3D_knee 

`thickness3D` run in knee mode (`-k`); there is no separate executable. As in
the 2D case, the only difference from the default cortex mode is where the two
boundaries come from: without `-k`, `thickness3D` derives them from the
white-matter/cortex tissue labels (`compute_boundary_cortex3D`); with `-k` it
instead reads a volume in which the cartilage boundaries are already encoded
(interior surface = 1, exterior = 0, cartilage = 128, rest = 255), as produced
by `compute_boundary_knee3D`. The Laplacian thickness computation itself is
identical.

```
Input: volume unsigned char file, where the interior boundary of the knee cartilage is 1, the exterior is 0, and the cartilage has label 128, and the rest has label 255  
Output: volume float file, representing the distances from one boundary to the other according to a metric defined by the laplacian field
```

# Examples
```sh
./thickness3D -n 10 -i 100 -k $DATA_DIR/knee-images/thickness/010/output.vol $DATA_DIR/knee-images/thickness/010/output_thickness3D.volf 512 512 70 
```
--------------------------------
# laplace2D 

Solves the 2D Laplace equation over a domain to produce the harmonic scalar
field that underlies the thickness algorithms: 0 on one boundary, 1 on the
other, smoothly interpolated in between. It is useful on its own for inspecting
or visualizing that field. The interior boundary is detected automatically from
the input domain.

## Usage
```
laplace2D [-c color] input.png output.png iterations [lambda]
              -c color output (0: gray, 1: red-blue, 2: random)
```
The domain should be a non-zero region (e.g. label 2) over a zero background,
so the interior boundary can be detected. The image dimensions are read from the
input PNG.

Output: filename ending in ".png" writes a PNG of the field (grayscale linearly
        normalized min..max to 0-255, or RGB when -c is 1 or 2); any other
        extension writes a raw float image, as before.

## Examples 2D:
```sh
./laplace2D data/domain_anillo_poisson.png output_laplace2D.png 100
./laplace2D data/domain_anillo_poisson.png output_laplace2D.flt 100
```
--------------------------------
# laplace3D 

Volumetric (3D) version of `laplace2D`: solves the Laplace equation over a
volume and writes the resulting harmonic scalar field as a float volume.

## Usage: 
```
laplace3D max1 max2 max3 input-prefix output3d-prefix iterations [lambda]
```

## Example 3D
```sh
./laplace3D 80 80 80 ../phantoms/thickness/domain_anillo_3d.vol output_laplace3D.volf 50
```
--------------------------------
# compute_boundary2D 

Extracts the inner and outer boundaries of a segmented 2D structure and writes a
domain image ready for the thickness tools. From a segmented image and its
original grayscale image, it isolates the structure of interest (selected by a
label and an intensity threshold) and labels the interior boundary as 1, the
exterior as 0, the interior region as 128 and the rest as 255. It is aimed at
knee-cartilage data and produces the kind of input `thickness2D_knee` expects.

## Usage
```
compute_boundary2D max1 max2 segmented_file.ush original_file.mri output2D.ush threshold
```

## Example: (knee data)
```sh
./compute_boundary2D 512 512 $DATA_DIR/knee-images/welsch/010/seg5/seg5-010.021 $DATA_DIR/knee-images/010/I_swap.021 ~/data/knee-images/thickness/010/output2D.chr 14 250 1
```

--------------------------------
# compute_boundary_knee3D 

3D boundary extractor for knee cartilage. From a segmented volume series and the
original grayscale series, it isolates the cartilage (by label and intensity
threshold) and produces the domain volume the thickness tools expect: interior
cartilage surface = 1, exterior = 0, interior = 128, outside = 255. This is the
volume you feed to `thickness3D -k`.

## Usage 
```
compute_boundary_knee3D [options] segmented_prefix original_prefix output3D.vols nrows ncols nslices
              -d (debug mode)
              -w (swapbytes) 
              -l label
              -t threshold 

   Input : segmented and original file series (unsigned short). We have to provide the label of the cartilage segmented to compute the boundaries on it.
   Output: volume unsigned char file, donde la superficie interior del cartilago tiene valor 1, la exterior tiene valor 0, lo de dentro tiene valor 128, y lo de fuera valor 255
```

## Example:
```sh
./compute_boundary_knee3D -l 14 -t 250 -w $DATA_DIR/knee-images/thickness/010/try5/input $DATA_DIR/knee-images/010/I ~/data/knee-images/thickness/010/output.vol 512 512 70
```

--------------------------------
# compute_boundary_cortex3D 

3D boundary extractor for the cerebral cortex. From a segmentation with a cortex
label and a white-matter label, it produces the domain volume with the inner
cortex surface = 1, the outer = 0, the interior = 128 and the rest = 255. It
plays the same role as `compute_boundary_knee3D` but is driven by tissue labels
instead of an intensity threshold. Note that `thickness3D` performs this same
boundary computation internally, so this standalone tool is mainly for
inspecting or reusing the boundary volume.

## Usage 
```
compute_boundary_cortex3D [options] input_prefix output3D.vols nrows ncols nslices
              -d (debug mode)
              -w (swapbytes) 
              -l label_wm 
              -c label_cortex 

   Input: serie de ficheros (unsigned short) que definen la segmentacion del cortex con etiqueta label_cortex, la segmentacion de la materia blanca con eqiueta label_wm, y el resto.   
   Output: Fichero de volumen (unsigned short) tipo vols, donde la superficie interior del cortex tiene valor 1, la exterior tiene valor 0, lo de dentro tiene valor 128, y lo de fuera valor 255
```

## Example
```sh
./compute_boundary_cortex3D -l 3 -c 2 -w $DATA_DIR/brain_124/Output $DATA_DIR/brain_124/Output_b.vol 256 256 124
```

--------------------------------
# Poisson2D

Solves the 2D Poisson equation over a domain (a variant of the Laplace solver
that includes a source term). This is experimental; the interesting output is
the map of local minima of the resulting field (`min_local.chr`).

## Usage
```sh
poisson2D max1 max2 input.chr output.pnm iterations [lambda]
```

## Example:
```sh
./poisson2D 256 256 domain_anillo_poisson.chr output_poisson.flt 100 
```

Lo interesante son los minimos locales (min_local.chr)
