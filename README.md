# mCSPC-NMA-publication
Code and data for mCSPC NMA publication

OpenBUGS must be installed to run these analyses. This is a separate software to R and RStudio and can be downloaded from the link below:
https://www.mrc-bsu.cam.ac.uk/software/bugs/openbugs/ 

# Instructions on running aggregate data models
These apply to all analyses i.e., base case proportional hazards, sensitivity analyses, class effects, inconsistency, safety, FACT-P and subgroup analyses. We focus on the OS base case analysis as an example.

Load the mCSPC_OS_doublet_ITT_base_case_FE script into your working space
Installation: Before using the R2OpenBugs package, you need to install it. You can install it from CRAN using the following R code:
install.packages("R2OpenBugs")

Before using R2OpenBugs, load the library in your R script:
library(R2OpenBugs)

If you haven't already installed the "readxl" package, be sure to install it from CRAN using the following command:
install.packages("readxl")

Load the readxl library in your R script: 
library(readxl)

Import the data from the excel file: Use the read_excel function to read the Excel file. Specify the file path as an argument. e.g., 
mCSPC_OS_doublet_ITT_base_case <- read_excel ("mCSPC OS doublet ITT base case.xlsx")

Running the Model:
Highlight the contents of the code and run the model.

Results:
After running the model, a HR cross effects csv file will be saved into your working directory. 

You can also access the results using the coda.samples function. This will provide a coda object with trace plots, summary statistics, and other useful information for posterior analysis.

# Alternative method to load aggregate data using the RStudio GUI
Alternatively, data can be easily loaded using the "import Excel data into R" option using the RStudio GUI
Click on "File" then choose "Import Dataset"
From the dropdown menu that appears, select "From Excel...". This will open a file dialogue where you can navigate to and select your Excel file.
Browse your computer and select the data file. Once you've selected the file, click "Open."

Specify Import Options:
RStudio will open a dialog box that allows you to specify import options. The dataset file contains one sheet and the first row is already set as column names so no need to specify any options.

Click the "Import" button. RStudio will read the Excel file and load it into a data frame in the R session.
