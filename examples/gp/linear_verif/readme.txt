1) Generate 60-run Latin hypercube design

	des_prep.r   : R code to generate design
	ranges       : Ranges of input parameters
	s-lhs.60.6   : Normalized design
	lhs.txt      : Design for input to GPMSA
        des_pre.pdf  : Bivariate projections of normalized design
        des_post.pdf : Bivariate projections of GPMSA design

2) Generate model output from lhs.txt

	y_mod = x1*beta1 + x2*beta2 + x3*beta3
        modout.r   : R code to generate model output y_mod
	y_mod.txt  : Model output to GPMSA

3) Construct experimental data matrix

	g.in       : Regression matrix
        y_1.dat    : Experimental data generated from g.in
	expdat.r   : R code to generate data input/output matrix
        y_exp.txt  : input/output matrix to GPMSA

4) Assumptions in GPMSA:

	a) yobs(ii).Sigy=0.05^2; (line 19 in sc.m)
	b) optParms.priors.lamOs.params=[64 4/25]; (line 16 in runmcmc.m)
	c) params.priors.theta.fname = 'gLogBetaPrior';
	   params.priors.theta.params = repmat([1 1],params.model.q,1);
           (lines 20-21 in runmcmc.m) 
        d) params.priors.lamVz.params=repmat([1 0.0001],params.model.lamVzGnum,1);
           (line 22 in runmcmc.m)
        e) params.priors.lamWs.params=repmat([1 0.0001],params.model.pu,1);
           (line 23 in runmcmc.m)

5) Output files:

	a) mcmc_scaled.txt   : All MCMC results, with calibration parameters on
				scaled domain ([0,1] for each input). Initial
				settings on first line of file.
	b) mcmc_unscaled.txt : All MCMC results, with calibration parameters on
				native domain. Initial settings on first line of file.

	Column order of parameters: theta1, theta2, theta3; betaV1, betaV2, betaV3;
                                    betaU1, betaU2, betaU3, betaU4, betaU5, betaU6;
                                    lamVz; lamUz; lamWs; lamOs; logLik; logPrior;
                                    logPost

6) Plots of results:

	verif1_cal_mcmc.dat : MCMC results from manual
	plotpost.r          : R code to plot marginal posteriors from GPMSA and manual	
        verif_1.pdf         : PDF plot of marginal posterior comparisons
