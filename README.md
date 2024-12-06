# mCSPC-NMA-publication
Code and data for mCSPC NMA publication

OpenBUGS must be installed to run these analyses. This is a separate software to R and RStudio and can be downloaded from the link below:
https://www.mrc-bsu.cam.ac.uk/software/bugs/openbugs/ 

# Instructions on running aggregate data models
These apply to base case proportional hazards, class effects, inconsistency, meta-regressions, safety, and subgroups analyses. We focus on the OS base case analysis as an example.

Load the mHSPC_OS_stampede_adjustment_FE script into your working space
Installation: Before using the R2OpenBugs package, you need to install it. You can install it from CRAN using the following R code:
install.packages("R2OpenBugs")

Before using R2OpenBugs, load the library in your R script:
library(R2OpenBugs)

If you haven't already installed the "readxl" package, be sure to install it from CRAN using the following command:
install.packages("readxl")

Load the readxl library in your R script: 
library(readxl)

Import the data from the excel file: Use the read_excel function to read the Excel file. Specify the file path as an argument.
mHSPC_OS_data_PO <- read_excel("~/Final NMA codes/Stampede adjusted/OS/mHSPC OS data PO.xlsx")

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

# Instructions on the fractional polynomial models

The fractional polynomials (PF) folder contains 2 subfolders: data and code.

The data folder contains the treatment codes (mHSPC OS trt codes.xlsx and mHSPC PFS trt codes.xlsx), which are only used for network plotting and to confirm that the correct treatment is used in the analyses. 

The FP_data_OS.csv and FP_data_PFS.csv files contain reconstructed individual patient data (IPD) from Kaplan-Meier (KM) curves, ready to be analysed. Each file is a data frame with 6 columns (s, a, dt r, z, time_j). s correspond to the study indicator, from 1 to 13. a is the correspondent arm in each study. dt is the time interval (in months) where a Piecewise Constant Hazard as approximate Likelihood is implemented. r is the number of events recorded during the time interval and z is the patients at risk at the beginning of the interval. time_j records the time at which the time interval ends (also in months).

The code folder implements the FP models on both PFS and OS outcomes. The OpenBUGS models are in the txt files, including fixed- and random-effects, first- and second-order FP models.

The OS_FP_models.R and the PFS_FP_models.R implement all the FP models explored during the analyses. OS_model_averaging.R and PFS_model_averaging.R implement model averaging from the 5 models with the lowest deviation information criteria (DIC) from the previous analyses. Finally, the OS_survival_curves.R and the PFS_survival_curves.R scripts implement survival extrapolation per treatment.

Please beware that the working directory has been handled using the `here` package. More details about this package can be found here: https://here.r-lib.org/
In order to make this package work first install the package with the install.packages("here") function. Then, load the package using the library(here) function. The package will detect the top-level of the project. To call any file use the here() function as follows: here("your folder name", "your file name"). 
