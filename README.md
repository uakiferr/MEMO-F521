# MEMO-F521
This study evaluated and compared DMR detection methods for a master's thesis in statistics.

The code provided here is based on the work of Jaffe et al. (2012), Butcher & Beck (2015), Mallik et al. (2019), Peters et al. (2020) and Alhassan et al. (2024). This study aimed to compare Bumhunter, Probe Lasso, DMRcate, and idDMR with simulated data and with real data based on the epi-signature found by Butcher et al. (2017). The simulated data was collected according to the Mallik et al. (2019) study and based on samples from the open-source Series GSE97362, accessible through the Gene Expression Omnibus link:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97362. 

The codes were modified slightly for this study. The computation and implementation were done on RStudio with R (v.4.3.3) on the Windows 11 Professional Operating System and a Dell Latitude 5590 (RP0FY) with 16 GB of RAM and Intel(R) Core(TM) i5-7300U CPU at 2.60GHz. So, all the parallel processing was deleted in order to run the codes. 

In this repository, the following files are available: 

- Simulation: All the R files are used to implement and compare the DMR detection methods.
- Real: All the R files used to implement and compare the detection methods of DMRs.
- Preprocessing: preprocessing of the datasets.
- Plots: All the R files to generate all the plots.
- Data: the type of data that was generated. 
