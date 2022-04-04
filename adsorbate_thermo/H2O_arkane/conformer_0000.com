%mem=5GB
%nprocshared=24
#P m062x/cc-pVTZ ! ASE formatted method and basis
opt=(calcfc,maxcycles=900,) freq IOP(7/33=1,2/16=3) scf=(maxcycle=900)

Gaussian input prepared by ASE

0 1
O                -0.0057000000        0.3852000000       -0.0000000000
H                -0.7961000000       -0.1947000000       -0.0000000000
H                 0.8018000000       -0.1905000000        0.0000000000


