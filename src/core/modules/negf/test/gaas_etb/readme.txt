Useful guide to set Negf calculations using Tight Binding and have less trouble

1) In the .geo file, make sure to define the quantum device length by integer multiples of the lattice constant lc.
    Example:
    | lc = 0.2504;
    | Npl = 6;
    | d = Npl*lc;
    | Point(1) = {0, 0, 0, 0.1}; # <--| This will be your quantum device region
    | Point(2) = {d, 0, 0, 0.1}; # <--|
    | [...]

2) In the `Quantum_contact` sections of the .tib file (inside the `Device` region), set the `length` option equal to lc defined in the .geo file.
    Example :
    | Device
    | {
    | [...]
    | Quantum_contact c1
    | {
    |     regions = < quantum_device_region_name >  
    |     length = 0.2504
    | }
    | 
    | Quantum_contact c2
    | {
    |     regions = < quantum_device_region_name > 
    |     length = 0.2504
    | }
    | [...]
    | }

3) In the `Solver` region of the `negf` module, set the `number_of_PL` option equal to Npl defined in the .geo file.
    Example:
    | Module negf
    | {
    | [...]
    | Solver
    | {
    |     [...]
    |     number_of_PL = 6
    |     [...]    
    | }
    | [...]
    | }
