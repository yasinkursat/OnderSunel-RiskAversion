# Risk-Aversion
"The role of risk aversion in a sovereign default model of polarization and political instability"
    
    1) Description of files for table and figure replications:
        a. Solving model for the benchmark parameterization with political risk (Column 3 of Table 2 in the manuscript): 
            i. Compile and run file political_risk_linear.f90. 
            ii. It takes 36 minutes and 47 seconds in Windows environment.
            iii. hpfilter.m: Matlab file that filters data using the Hodrick-Prescott filter.
            iv. simulate.m: Matlab file that produces the moments that are reported in the tables of the paper.
        b. Solving model for alternative parameterizations and to replicate the rest the columns in Table 2 as well as the moments in Table 3
            i.  Replace the corresponding parameters in the FORTRAN file. NOTE: Risk aversion parameter in the code is called "sigma" while in the paper it is represented with "gamma". 
        c. Solving the model for Online Appendix:
            i. To solve the discretized state space method; run political_risk_dss.f90 (first two columns of Table 1 of Online Appendix). Corresponding matlab file for the simulation is simulate_dss.m.
            ii. To solve the spline method; run political_risk_spline.f90 (last two columns of Table 1 of Online Appendix). Corresponding matlab file for the simulation is simulate_spline.m.
            iii. To solve for the long-term debt model; run political_risk_longterm.f90. Corresponding matlab file for the simulation is simulate_longterm.m

    2) How to run the codes

	FORTRAN files require having an access to the IMSL library from which the files invoke some subroutines from. Save your files in a subdirectory named "graphs" for each FORTRAN run. 
	Save matlab files in the same subdirectory. Necessary seeds for random number generation are all specified in the corresponding FORTRAN files. 
	We use a random number generator to draw sequences of income shocks, probability of exclusion and probability of office change. 
	We keep the draws so that the same values can be used for each sample. Routine RNSET of IMSL is used to initialize the seed, which is set to be 139719.
