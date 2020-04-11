# pyroproQ3dock

                                ____   _____     __           __  
   _________  ____  _________  / __ \ |__  /____/ /___  _____/ /__
  / ___/ __ \/ __ \/ ___/ __ \/ / / /  /_ </ __  / __ \/ ___/ //_/
 / /  / /_/ / /_/ / /  / /_/ / /_/ / ___/ / /_/ / /_/ / /__/ ,<   
/_/   \____/ .___/_/   \____/\___\_\/____/\__,_/\____/\___/_/|_|  
          /_/                                                     

RoproQ3dock v 0.2

1. This is a program which takes SabDab chothia formatted crystal structures of protein antigen complexes. (Single or in a Folder)
2. It cleans the structure choosing the best complex by Bfactor and relabelling and renumbering the chains (L,H,A).
3. Single chain antibodies and more than three chain systems are eliminated.
4. The complex is minimized, packed with pyrosetta and then randomized.
5. Rigid body FFT docking is performed on the GPU.
6. The top 96 scored structures are further relaxed and interface energy is analyzed.
7. ProQ3D scoring is perfomed on the complexes
8. Scorefile is returned ranked by rmsd to input reference
