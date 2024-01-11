The files included here include source code of the Bryan Cloud Model CM1, release 19.5.  The individual original modules can be downloaded here: https://www2.mmm.ucar.edu/people/bryan/cm1/.
The present version includes modifications mainly in the parcel routine (parcel.F), where diagnostics have been added that calculate the rhs of the vorticity equation along the parcel
trajectories.  

Files included:

1. onefile.F: All relevant, precompiled modules linked together
2. parcel.F: modified parcel-trajectory code that includes vorticity diagnostics
3. namelist.input: Namelist settings to simulate a dry CBL with dust-devil-like vortices
4. input_sounding: Slightly stable, dry base-state sounding  
