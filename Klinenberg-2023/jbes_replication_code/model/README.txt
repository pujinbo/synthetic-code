The model is composed of two:

bitto_model_revised.R: BL-TVP, the proposed model, utilizing the MCMC algorithm proposed in Bitto-Fuhrwerth Schnattner (2019)

bscm_horseshoe.stan: The .stan code for BSCM-Horseshoe, one of the models BL-TVP is compared to. This is used in the TVP simulation and main results. The time invariant simulation uses the `bayesreg` R package. After conferring with the authors of BSCM-Horseshoe, it became clear both programs perform the same. `bayesreg` is sufficiently faster.