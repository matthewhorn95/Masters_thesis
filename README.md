# Masters_thesis

This is some work from my master's thesis which involved phenological modeling of Finnish moths applying methods from survival analysis. I am mainly putting this here as an example of my abilities and will include the pdf of my master's thesis as well as the grade I received. 

I didn't save some of the files, so this is not the entirety of the project, and the modeling file is not the final one, since the most recent version I could find was from March 1st. Still, it should show some of what I'm capable of, even if some of the exploratory analysis and data processing is missing.

The code will not be able to be run, because I cannot share the data, as it is the property of the Finnish government and is confidential.

There is a R markdown file in which I extracted some of the covariates at the sampling sites, an R Markdown file in which I did some data processing, EDA and the modeling and predictions, and a STAN file which controls the actual model itself, intaking a data vector created in the modeling file and computing the joint posterior distribution as MCMC samples using No U Turn Sampling. This is all detailed in the thesis pdf itself.
