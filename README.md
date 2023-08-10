# melt
set of lammps input files and python analysis scripts to calculate the melting point of 70-Fe 21-Ni 9-Cr

the melting point is calculated by placing a solid phase and a liquid phase next to each other. below the melting point, the solid phase will survive. above the melting point, the liquid phase will survive instead. we can detect the melting point by running simulations at some guessed temperatures, and seeing at which point the behavior changes.

to run the LAMMPS input at a given control pressure P (in bars) and a temperature T (in kelvin), run:

``lmp -var temperature T -var pressure P -in melt.in``

replace `lmp` with whatever your LAMMPS executable is named. for example, if I want a run at 1 bar and 1000 K and my LAMMPS executable is named `lmp_win.exe` I would run:

``lmp_win.exe -var temperature 1000 -var pressure 1 -in melt.in``

of course, you're free to use any other flags, like `-sf gpu` if you're using a GPU-built version of LAMMPS

this simulation will have some output files. the important output file is named `equil_T_P.dump`, which has the phase-phase competition dynamics stored. be careful with renaming these files - the analysis script assumes the files are named this way.

then, to run the analysis script, run:

``python percent_solid.py T_low T_high T_step P``

here `T_low` is the lower temperature (in K), `T_high` is the upper temperature (in K), `T_step` is the step between these two temperatures (in K), and `P` is the control pressure (in bars)

this script will generate a plot of percent fcc vs. temperature and save it to `melting_P.png`. if the run was long enough, there will be a very clear transition near the melting point which one can visually identify

this python script depends on the external packages ovito and matplotlib. both of these can be installed with `pip`

uncomment line #7 if you get a weird error with the matplotlib import
