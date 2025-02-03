# project3
Start date: 10/29/2024
end date: 02/03/2025
The problem need to be fixed:
1. the scf energy result is different from the reference. What I got is  "The initial SCF electronic energy is:  -129.84540364532984", however the reference is "The initial Hartree-Fock electronic energy for the H2O test case is -125.842077437699 Hartrees." Even though the final SCF energy is only differed at the level of 0.1. But the result from pyscf is differed at the level of 0.01. This need to be furthur investigated.
2. I used a lot global variables, not sure whether it is a good idea or not. I will try to fix this problem in the future.

     