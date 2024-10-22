# Coding for Mordern Quantum Chemistry.
Current goal 
- Understand Restricted Hartree Fock Calculations better with the code ready to use.
- Extend the current code to Finite-Temperature Restricted Hartree Fock Calculations
- Extend the basis set constructions "STO-3G, STO-6G, ..." <br>
Edited date: 04022024 <br>
Hopefully, this can be finished by the start of the summer:)

Updated date: 06122024
Trying to get a module for gaussian function, where I can normalize the function, and calculate all the one- and two-electron integral.

Update on 10162024: nothing have been done:(

Update on 10222024: 
But I will try to follow the course: https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2301
And project 1 have been finished. The archtecture used is fortran package manager(fpm)
'fpm test' can get the similar results as shown in the course website. 
Question to be further investigated:
1. what is the difference between subroutine and function? 
    I know the fact that funciton result in one values and subroutine result in multiple modified values, but how to use them? Or, more specifically, how to call a subroutine from another subroutine properly and how to share the data globally and locally?
2. printing function is poor coded.
