README FILE FOR NUMERICAL SECTIONS OF <PAPER NAME>
==================================================

*Simulations*
lmselect_gamma_sim_marginal_mcore.R: linear regression simulation, generalized gamma bootstrap
	outputs stored in lmOutMarginalGamma.Rda
	plotting done using Sim_plotsMarginalGamma.R
	plots saved in simplot_gammax.pdf, x=2,4,6,8	
lmselect_moon_sim_marginal_mcore.R: linear regression simulation, moon bootstrap
	outputs stored in lmOutMarginal.Rda
	plotting done using Sim_plotsMarginal.R
	plots saved in simplotx.pdf, x=2,4,6,8
lmselect_wild_sim_marginal_mcore.R: linear regression simulation, wild bootstrap
	outputs stored in lmOutMarginalMooN1.Rda
	plotting done using Sim_plotsMarginalMooN.R
	plots saved in simplotmoonx.pdf, x=2,4,6,8
lmmselect_wild_sim_marginal_gamma.R: linear mixed model simulation, generalized gamma bootstrap
	outputs stored in lmmout150marginalgamma.Rda, lmmout600marginalgamma.Rda
lmmselect_wild_sim_marginal_mcore.R: linear mixed model simulation, wild bootstrap (not reported in paper)
	outputs stored in lmmout150marginal.Rda, lmmout600marginal.Rda


*Indian Monsoon data analysis*
monsoon_lmmselect_gamma_parallel.R: performs model selection on monsoon data, produces bestmodel.Rda which is not supplied because of large size (161 mb)
monsoon_lmmselect_gamma_plots.R: generates the 4 plots rolling_xyz.pdf from bestmodel.Rda

