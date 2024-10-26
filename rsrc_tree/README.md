## Coding Principles

1. Every file in the orig `cvxpy` source tree must have a corresponding file,
   even if empty file. 3. The first line should always be the source, e.g.: file
   `atoms/von_neumann_entr.py` will have `## CVXPY SOURCE:
   cvxpy/atoms/von_neumann_entr.py` as the first line.
2. Newly created files need to be added to the `rsrc_tree` (see below
   also) in the appropriate place. No files in the `R` directory
   should be edited and indeed, all files there are ephemeral and
   generated from the `rsrc_tree`. 
2. R-specific code sections should ideally be delineated using a block
   thus:
   ```
   ## Begin R-specific Code Section
   ...
   ## Begin R-specific Code Section
   ```
3. Auxillary R code needed to make things work go in `zzz_R_specific`
   directory with appropriate comments. Cross links should be made as
   needed to enable anyone to follow the code in its entirity.
4. All code must be `lint`ed using `lintr` before PR so that we don't
   have any silly diffs showing up merely because of typographic
   (extra spaces etc.) or stylistic differences. 

## Package Building Principles
1. All files used for R code will reside in `rsrc_tree` mirroring the
   structure of `cvxpy/cvxpy`.
2. The `R` directory for the package is generated using a script
   `inst/copy_source.R`. 
3. The process for documenting the package first runs the script in
   step 2 and then executes `devtools::document()`. 
   
   
