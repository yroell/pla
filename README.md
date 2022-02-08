# Parallel Line Assay (PLA)
The code for conducting a PLA analysis has been provided with a visualization of the major steps within the code. Below talks about how to implement the PLA code and discusses the steps occurring within the script.

### Implementing Code:
The code has been designed so only lines 13-16 have to be changed. Once those lines have been changed to point to the specific data on your computer, highlight all of the code (Windows: Ctrl + A; Mac: Command + A) and click Run (alternative – Windows: Alt + Enter; Mac: Command + Enter). The folder structure will be created on your computer in the path indicated (line 13) and the data will be read (using the name indicated on line 14; needs to be a csv file). The folders and output files will use the lab name (line 15) within the naming convention. The standard ID (line 16) lets the script know what everything should be compared to (as the reference).

### Adjustments to your dataset
The csv columns need to have columns with these names: **Result**, **Dilution_Factor**, **Sample_ID**, **Antigen**, and **Antibody**. The standard ID set on line 16 needs to be within the **Sample_ID** column.
- If the data trying to be read in contains nuances not incorporated within the script pertaining to the naming conventions within the antibody, antigen or result column, the script might crash. 
- The best practice would be for the antibody and antigen columns to only contain a single word and no special characters. If multiple words are needed, separate words with an underscore. 
   - For example, instead of “Total (IgM/IgG/IgA)”, “Total” would be ideal, and instead of using a percent sign (%) with neutralization, “percent_neutralization” would be ideal. The last bit of code within this section creates the output logfile to begin keeping track of what steps are occurring for each sample.

### Steps in Code:
Once the user inputs have been specified, the rest of the code will run without anything needing to be changed. This section of the text will breakdown the code into smaller pieces to explain what is occurring.

*Lines 20-34*
Packages used within the script (stringr and dplyr) are installed and/or loaded. The next chunk of code is a function to clean up how numbers on the figures are presented.

*Lines 36-60*
The path set in line 13 is the <ins>**working directory**</ins>. A folder is created that will have the output files with the results located here when the script is complete and contains a subfolder where the figures will be. This figure folder will also have a subfolder where the non-parallel samples will be if it is possible to plot them as well. 
The output files are: 
1. A csv that contains all the unique samples (sample ID, antigen, and antibody combination) with the relative slope, unknown BAU, and confidence interval.
2. A text file that documents every step that is occurring to help with debugging if the code does not finish (every step within the code is recorded here).
3. A text file that contains which samples were not plotted and the reason why. The csv file with the raw data is read in next and some cleaning of the data occurs.

*Lines 62-89*
The <ins>**dilution factor**</ins> column is used to determine what is the correct format to use as the dose. The <ins>**dilution factor**</ins> can be a character that has a numerator, forward slash, and a denominator. For example, “1/50” or “1/300” so the dose will be 1 divided by 50 or 1 divided by 300. The dilution factor can be an integer (whole number above 1) which means the dose will be 1 divided by the integer provided. The <ins>**dilution factor**</ins> can also a float (decimal between 0 and 1) which means the dose will be the float. The next step log transforms the dose and result column. The last bit of code in this section finds all of the unique combinations within the different <ins>**sample ID**</ins> (excluding the <ins>**standard ID**</ins>) and antigen-antibody combinations.

#### *Lines 91-553*
This is the bulk of the code and will be examined in smaller chunks of code but this all repeats for each unique combination identified in the lines of code beforehand. The next portions explained will assume we are only working with one unique combination.

*Lines 93-143*
The data is split so that we have a dataset that now only contains the samples relating to the <ins>**standard ID**</ins>. Data is removed if there are not three replicates from the same dilution. If a <ins>**standard ID**</ins> and antigen-antibody combination does not exist, this is recorded in the logfile and the script moves on to the next unique combination. If the combination exists, a line is created from the <ins>**standard ID**</ins> samples and the slope and intercept of the line are obtained.

*Lines 145-215*
The top and bottom dilutions are removed as necessary. If the top or bottom 3 points are not within 20% of the previous slope obtained and there are 5 or more dilutions and the R2 of the line created is less than 0.95, the top or bottom point is removed. This is repeated until either the slope is within 20% or there are less than 5 dilutions or the R2 is 0.95 or greater. Every time a point is removed, a new line is created. 
- *The removing of points is to clean up the tail ends from a sinusoidal curve and trying to get the straight part from the curve.*

*Lines 217-293*
The data is now split so that we have a dataset that only contains the samples with the unique combination identified previously and points are removed if they fall outside of the top or bottom dilutions within the standard dataset or if they have less than three replicates. If this <ins>**sample ID**</ins> and antigen-antibody combination does not exist, this is recorded in the logfile and the script moves on to the next unique combination. If the dataset contains less than three dilutions, this is recorded in the logfile and the script moves on to the next unique combination. If the combination exists, a line is created from this dataset and the slope and intercept of the line are obtained.

##### Script diverges depending on the slope of the lines
The next few chunks of codes either plot the data if they are already parallel or tries to make the data parallel as best as possible. If the data cannot be parallel after all the steps, the script moves on to the next unique combination.

*Lines 295-373: data is already parallel*
If the <ins>**standard ID**</ins> dataset and the <ins>**sample ID**</ins> unique combination is already parallel, meaning the <ins>relative slope is within 20% of each other (relative slope = standard ID slope / unique combination slope)</ins>, no cleaning is required and the BAU can be calculated with the data plotted. Bootstrapping is implemented to calculate the BAU and confidence interval. The relative potency is calculated 100 times (based on 100 datasets created by sampling with replacement at each dilution to get 100 different lines for both the standard ID and unique combination datasets). 
The potency is calculated by: 
1. Taking 100 points between the minimum and maximum result from the standard ID dataset. 
2. Subtracting the standard ID intercept and dividing by the standard ID slope to get the standard ID potency (repeated for the unique combination as well to get the unique combination potency).
3. The log relative potency is the standard ID potency minus the unique combination potency. The log relative potency is back transformed and averaged over the 100 points and multiplied by 1000 to get the relative potency. 
- The BAU is the average of the 100 replicas of the relative potencies. The BAU and confidence interval are recorded and a plot is created.

*Lines 376-544: data is not parallel*
Some of the <ins>**standard ID**</ins> and unique combination datasets are not parallel (relative slope is greater than 20%), so cleaning is required. The next chunk of code tries to make the datasets parallel as best as it can. This is done by removing either the top or bottom dilution and calculating the relative slope. Removing whichever point that creates a new line that is closer to or less than 20% for the relative slope is the point removed. This is repeated until the relative slope is less than 20% or there are fewer than three dilutions left. If there are fewer than three dilutions remaining, it is noted that the data cannot be parallel so the script moves on to the next unique combination. If removing points ever makes the relative slope less than 20%, the BAU is calculated and the data is plotted just as it was done in lines 295-373.

*Lines 547-549: data cannot be made parallel*
The unique combination cannot be made parallel to the <ins>**standard ID**</ins> after all the steps taken in lines 376-544. This is recorded and the script continues to the next unique combination.

*Lines 544-559*
The BAU and confidence intervals recorded are saved and the logfile is saved as well.

*Lines 563-585: commented out by default*
If running this script from multiple datasets (data from different labs), this chunk of code will read in all of the data created from running this code and create a csv file that contains all of the different datasets into one large dataset. Only need to uncomment and run if one large dataset is needed (i.e. further analysis required across the labs).
