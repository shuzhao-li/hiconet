
"""
BTM is now on pypi.

Example to convert gene level data to BTM module activity scores:
>>> genetable_to_activityscores('MCV4_D3v0_genes.txt', 'MCV4_D3v0_BTMactivity.txt')

The output file is MCV4_D3v0_BTMactivity.txt. The module activity scores are computed as the mean
value of member genes. You can use these activity scores to perform further statistical test of your
choice.
"""

from btm.btm_tool import genetable_to_activityscores


if __name__ == '__main__':
    import sys
    infile = sys.argv[1]
    outfile = 'BTMactivity_' + infile
    genetable_to_activityscores(infile, outfile)
    print("Result written to " + outfile)
