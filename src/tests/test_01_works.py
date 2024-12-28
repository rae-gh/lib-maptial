import os, sys
from pathlib import Path
sys.path.append(os.path.join(os.path.dirname(Path(__file__).parent)))


DATADIR = "src/data/"
ls_structures = ['6eex']

def test_01():    
    print("test_01 backbone")
    

if __name__ == "__main__":    
    test_01()