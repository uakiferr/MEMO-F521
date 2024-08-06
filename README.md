# MEMO-F521
Evaluation and comparison of DMRs detection methods for a master's thesis in statistic

The code provided here is based on the work of Jaffe et al. (2012), Butcher & Beck (2015), Mallik et al. (2019), Peters et al. (2020) and Alhassan et al. (2024). The aim of this study was to compare Bumhunter, Probe Lasso, DMRcate and idDMR with simulated data and with real data based on the epi-signiture found by Butcher et al. (2017). The simulated data was done according to Mallik et al. (2019) study and based on samples from the open source Series GSE97362, accessible though the Gene Expression Omnibus link:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97362. 

The codes were a bit modified for this study. The computation and implementation was done on RStudio with R (v.4.3.3) on the Windows 11 Professional Operating System and a Dell Latitude 5590 (RP0FY) with 16 GB of RAM and Intel(R) Core(TM) i5-7300U CPU at 2.60GHz. So all the parrallel processing was delete in order to run the codes. 

In this repository, you can find the following files: 

- Simulation : All the R files used to implement and compare the DMRs detection methods.
- Real : All the R files used to implement and compare the DMRs detection methods.
- Preprocessing : preprocessing of the datasets.
- Plots : All the R files to generates all the plots.
- Data : the type of data that I generated. 
