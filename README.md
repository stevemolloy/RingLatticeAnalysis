# Analysis of an accelerator lattice

## Quick start
```console
$ make
$ ./bin/rla ./lattices/m4U_240521_b03_03_07_06.mad8 -E 3.0 -h 176 -p 20 --save_disp dispersion_output.csv
```
# Usage
The `rla` executable is a command-line application that takes a MAD8 lattice file as input and generates some analysis outputs of the lattice.

If the analysis detects that the lattice is closed (that is, the total bending angle of the line multiplied by the periodicity is very close to 360 degrees)
it will perform some extra analysis based on periodicity conditions.

## Flags
- `-E <total energy>` -- the total energy of the particles
- `-h <harmonic number>` -- the harmonic number of the ring. Used to calculate the RF frequency.
- `-p <periodicity>` -- the periodicity of the line. That is, the number of times the line repeats to make the total line.
- `--save_disp <filename>` -- save the calculated dispersion values to `filename`.
