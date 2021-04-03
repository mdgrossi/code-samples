# code-samples
Sample R and Python scripts to demonstrate data processing and analysis

These scripts are intended to be short samples to demonstrate my capabilities as a scientist at handling, processing, and analyzing large data sets using various statistical methods pertinent to different research projects I have conducted over the years. Each script represents a glimpse into much larger routines and/or pipelines of routines that collectively carry out different research projects in oceanography.

SYNOPSIS

**qc_summary.R** is one component from a pipeline of routines that process and analyze global Argo data. This script compiles QC flags while the next script in the routine applies these flags to filter out the undesired data. This was created to handle >2 million Argo profiles (from the program's conception through 2015). For more information, see the docstring in the file.

**autoarima.R** was designed to automatically and efficiently fit ARIMA (AutoRegressive Integrated Moving Average) models to thousands of observed ocean Lagrangian drifter velocity time series. The **R** function *auto.arima()* eliminates the need to conduct a manual parameter search for all time series and implementing it with **R**'s *apply* family of functions allows this to be done in a vectorized fashion. The parent script (not provided) compiles many such forecasts from the defined function *run.autoarima()* and writes them out to a netCDF file.

**ssa.py** demonstrates an object-oriented approach to implementing univariate and multivariate singular spectrum analysis (SSA) in python. As noted and credited in the docstring, this started with a simple online demonstration that was motified and expanded to suit the needs of the project at hand. The working jupyter notebook shows parts of the script in use.

The annotations contained in these scripts are standard practice in my coding career. It helps me review and recall what I wrote when looking back after some period of time, but it is not just for me. I am a team player who values reproducability and sharing science; as such, I believe any code I write for any project should be interpretable by anyone who needs to look at it. While this is usually common practice in computer science, it is not so in physical sciences.

I am proficient in both **R** and Python and generally select the language based on the task at hand: being a statistical programming language, **R** is usually best for statistical analysis and visualization and often has nifty tools that haven't been developed (or at least not publically) in other languages, such as *auto.arima()*. But while machine learning (ML) can be done in **R** as well as Matlab and many other languages, the ML community has largely settled on Python for its readability, simplicity, and open source status. 


CONCLUDING REMARKS

New scripts or notebooks may be added as needed or appropriate. Interests are encouranged to reach out to me with questions, requests, or solitications.

For an example of a "just for fun" side project I took up, see the respository **covid19**, which contains plots of COVID-19 data from various parts of the country (regions were selected based on the location of family and friends.) These plots are/were automatically updated daily throughout the global coronovirus pandemic.


LEGALITY NOTICE

The contents of this respository are provided demonstration purposes only. All rights to the products herein are reserved under standard United States copyright law. No source codes may reproduced, distributed, or used to create derivative works without expressed permission by the owner of this respository.
