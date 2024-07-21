## Principles

1. Every file in the orig source tree must have a corresponding file,
   even if empty file.
2. The first line should always be the source, e.g.: file
   `atoms/von_neumann_entr.py` will have `## CVXPY SOURCE:
   cvxpy/atoms/von_neumann_entr.py` as the first line.
3. Auxillary R code needed to make things work go in `zzz_R_specific`
   directory with appropriate comments. Cross links should be made as
   needed to enable anyone to follow the code in its entirity.
   
   
