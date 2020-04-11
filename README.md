
    ________________________________________________________________________________________________________
    ___________________________________________________________/ __ \__|_   /______/ /_________________/ /__
    ___/____\__/ /_/ /__/ ___/_/ __ \___/ __ \__/ ___/_/ __ \_/ /_/ /___/  <______/ /__/ __ \_/ ___/__/ //_/
    __/ / / /_/ /_/ /__/ /____/ /_/ /__/ /_/ /_/ /____/ /_/ // /_/ /_____/ /_/ /_/ /__/ /_/ // /_____/   <__
    _/ .___/__\__, /__/_/_____\____/__/ .___/_/_/_____\____/_\___\_\_/____/__\__,_/___\____/_\___/__/_/|_|__
    / /______/____/__________________/_/____________________________________________________________________

# pyroproQ3dock

RoproQ3dock v 0.2

1. This is a program which takes SabDab chothia formatted crystal structures of protein antigen complexes. (Single or in a Folder)
2. It cleans the structure choosing the best complex by Bfactor and relabelling and renumbering the chains (L,H,A).
3. Single chain antibodies and more than three chain systems are eliminated.
4. The complex is minimized, packed with pyrosetta and then randomized.
5. Rigid body FFT docking is performed on the GPU.
6. The top 96 scored structures are further relaxed and interface energy is analyzed.
7. ProQ3D scoring is perfomed on the complexes
8. Scorefile is returned ranked by rmsd to input reference
