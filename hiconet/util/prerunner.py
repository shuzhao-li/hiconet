"""prerunner.py

Check data matrix for shifted 1st row as in R convention.
-- This is potentially a bad idea. Better to check at i/o when reading file?









To-do:

Check missing values

Check annotation mapping



➜  test python ~/Code/hiconet/prerunner.py SDY80_pH1N1_2009_Transcriptomics.txt 
Header shifted and a new file written.
"""

import os

def check_and_fix(f, workdir='.', delimiter='\t'):
    s = open(f).read()
    [line1, line2] = s.splitlines()[:2]
    N1, N2 = len(line1.split(delimiter)), len(line2.split(delimiter))
    if N1 == N2:
        print("Header length = data column. No fix needed.")
    elif N1 + 1 == N2:
        with open(os.path.join(workdir, 'shiftcorrected_'+f), 'w') as O:
            O.write('cell1\t' + s)
        print("Header shifted and a new file written.")
    else:
        print("Mismatched length? No action taken.")





if __name__ == '__main__':
    import sys
    check_and_fix(sys.argv[1])
