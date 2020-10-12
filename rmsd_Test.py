import os
import re
import sys
import fileinput
import subprocess
from pathlib import Path
import docking_structure_randomizer
import anarci_to_pdb
import parapred_to_blocking
import decoygen
import zrank_processing
import rmsd
import numpy as np
import pandas as pd


test = rmsd.main('6BGT_1.pdb','complex_156.pdb')
print(test)