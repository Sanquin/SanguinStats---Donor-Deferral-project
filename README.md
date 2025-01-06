# R script for alternative deferral

Generic R-code was developed that allows performing an initial analyses of potential benefits of a mean Hb-level based donor deferral policy that compensates for biological and measurement variability. A more complete discourse of this approach is provided in the paper *“Why the Majority of On-Site Repeat Donor Deferrals Are Completely Unwarranted…”* (Transfusion 2022, 62 (10), 2068–2075. https://doi.org/10.1111/trf.17085).

## Packages
Need `dplyr` and `zoo` and `caret` for the alternative deferral analysis.

## Data
The input required to run the code is an .rds file containing 4 variables:

| Variable name | Variable description | Variable data type |
| --- | --- | --- |
| KeyID | Unique identifier for each donor | integer | 
| Sex | indicator for donor being male (M) or female (F) | character |
| DonDate | date of donation | date | 
| Hb | donor Hb at donation | numeric | 

If sufficient data is available it is recommended to include only donors for which the full donation history is available. For blood establishments with pre-donation and post-donation screening, we recommend keeping ONLY the post-donation screening result. 

## Parameters in the analysis code

In the codefile (*alt_def.R*) a number of parameters need to be specified by the user. Each of these parameters are to be stored in a variable:
1)	*FILE_DIR*: the directory in which the .rds file is located
2)  *DATAFILE_NAME*: the name of the .rds file
3)	*CUTOFF_M/CUTOFF_F*: Minimum acceptable Hb levels for males and females in the units of measurement also used in the data file (e.g. g/dL, g/L, mmol/L, etc.)
4)  *UNITS*: indicates the unit of measurement in the data file as a string (e.g. g/dL, g/L, mmol/L)
4)	*TEST*: indicates whether a test version of the code is ran (using only 10,000 records). Set this to F for the full analysis.


