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


Alhassan, D., Olbricht, G. R., & Adekpedjou, A. (2024). Differential methylation region detection via an array-adaptive normalized kernel-weighted model. PLOS ONE, 19(6), e0306036. https://doi.org/10.1371/journal.pone.0306036
Butcher, D. T., Cytrynbaum, C., Turinsky, A. L., Siu, M. T., Inbar-Feigenberg, M., Mendoza-Londono, R., Chitayat, D., Walker, S., Machado, J., Caluseriu, O., Dupuis, L., Grafodatskaya, D., Reardon, W., Gilbert-Dussardier, B., Verloes, A., Bilan, F., Milunsky, J. M., Basran, R., Papsin, B., … Weksberg, R. (2017). CHARGE and Kabuki Syndromes: Gene-Specific DNA Methylation Signatures Identify Epigenetic Mechanisms Linking These Clinically Overlapping Conditions. The American Journal of Human Genetics, 100(5), 773–788. https://doi.org/10.1016/j.ajhg.2017.04.004
Butcher, L. M., & Beck, S. (2015). Probe Lasso: A novel method to rope in differentially methylated regions with 450K DNA methylation data. (Epi)Genomics Approaches and Their Applications, 72, 21–28. https://doi.org/10.1016/j.ymeth.2014.10.036
Jaffe, A. E., Murakami, P., Lee, H., Leek, J. T., Fallin, M. D., Feinberg, A. P., & Irizarry, R. A. (2012). Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies. International Journal of Epidemiology, 41(1), 200–209. https://doi.org/10.1093/ije/dyr238
Mallik, S., Odom, G. J., Gao, Z., Gomez, L., Chen, X., & Wang, L. (2019). An evaluation of supervised methods for identifying differentially methylated regions in Illumina methylation arrays. Briefings in Bioinformatics, 20(6), 2224–2235. https://doi.org/10.1093/bib/bby085
Peters, T. J., Buckley, M. J., Statham, A. L., Pidsley, R., Samaras, K., V Lord, R., Clark, S. J., & Molloy, P. L. (2015). De novo identification of differentially methylated regions in the human genome. Epigenetics & Chromatin, 8(1), 6. https://doi.org/10.1186/1756-8935-8-6
Peters, T. J., Buckley, M. J., Yunshun Chen, Smyth, G. K., Goodnow, C. C., & Clark, S. J. (2021). Calling differentially methylated regions from whole genome bisulphite sequencing with DMRcate. Nucleic Acids Research, 49(19). https://academic.oup.com/nar/article/49/19/e109/6329576
