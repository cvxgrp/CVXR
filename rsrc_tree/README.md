## Coding Principles

1. Every file in the orig `cvxpy` source tree must have a corresponding file, even if empty file. The first line should always be the source, e.g.: file `atoms/von_neumann_entr.py` will have `## CVXPY SOURCE: cvxpy/atoms/von_neumann_entr.py` as the first line.
2. Newly created files need to be added to the `rsrc_tree` (see below also) in the appropriate place. No files in the `R` directory should be edited and indeed, all files there are ephemeral and generated from the `rsrc_tree`. 
3. R-specific code sections should ideally be delineated using a block thus:
   ```
   ## Begin R-specific Code Section
   ...
   ## Begin R-specific Code Section
   ```
4. Auxillary R code needed to make things work go in `zzz_R_specific` directory with appropriate comments. Cross links should be made as needed to enable anyone to follow the code in its entirity.
5. All code must be `lint`ed using `lintr` before PR so that we don't have any silly diffs showing up merely because of typographic (extra spaces etc.) or stylistic differences. We may need to enforce this via Github Actions for our sanity
6. 

## Package Building Principles

1. All files used for R code should reside in `rsrc_tree` mirroring the structure of `cvxpy/cvxpy`.
2. The collation order of files in the `rsrc_tree` into the generated `R` folder matters for the R package to be build-able. This order is dictated in the file `inst/all_files.csv` and is used by the script `inst/copy_source.R` to populate a new `R` folder whenever it is run. So any changes made to files in the `R` directory are ephemeral and will be ruthlessly overwritten.
3. 
3. To document the package, first run the script in step 2 and then execute `devtools::document()`. 

   
   
