
                                                                ____   ______       __                  __
    ___________________________________________________________/ __ \__|__  /______/ /_________________/ /__
    ___/ __ \__/ /_/ /__/ ___/_/ __ \___/ __ \__/ ___/_/ __ \_/ / / /___/_ <______/ /__/ __ \_/ ___/__/ //_/
    __/ /_/ /_/ /_/ /__/ /____/ /_/ /__/ /_/ /_/ /____/ /_/ // /_/ /_____/ /_/ /_/ /__/ /_/ // /_____/   <__
    _/ .___/__\__, /__/_/_____\____/__/ .___/_/_/_____\____/_\___\_\_/____/__\__,_/___\____/_\___/__/_/|_|__
    /_/      /____/                  /_/

# pyroproQ3dock

PyRoproQ3dock v 0.2

1. This is a program which takes SabDab chothia formatted crystal structures of protein antigen complexes. (Single or in a Folder)
2. It cleans the structure choosing the best complex by Bfactor and relabelling and renumbering the chains (L,H,A).
3. Single chain antibodies and more than three chain systems are eliminated.
4. The complex is minimized, packed with pyrosetta and then randomized.
5. Rigid body FFT docking is performed on the GPU.
6. The top 96 scored structures are further relaxed and interface energy is analyzed.
7. ProQ3D scoring is perfomed on the complexes
8. Scorefile is returned ranked by rmsd to input reference

# Requirements:

1. pyrosetta (https://www.rosettacommons.org/software/license-and-download)
2. proQ3 (https://bitbucket.org/ElofssonLab/proq3/src/master/)
3. megadock (https://www.bi.cs.titech.ac.jp/megadock/)
4. pdb-tools (https://pypi.org/project/pdb-tools/) 
    currently these need to be in the PATH
5. Parapred  (https://github.com/eliberis/parapred)
    also needs to be in the PATH
6. Anarci (http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/)
    also needs to be in the PATH
7. Packages - All the packages required to run the above programs plus  
    a. biopython  
    b. numpy  
    c. pandas  
    d. biopandas  
    e. mpi4py  
    f. GNU parallel
