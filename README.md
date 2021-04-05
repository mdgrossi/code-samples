# code-samples
Sample R and Python scripts to demonstrate data processing and analysis

These scripts are intended to be short samples that demonstrate my capabilities as a scientist at handling, processing, and analyzing large data sets using various statistical methods pertinent to different research projects I have conducted over the years. Each script represents a glimpse into much larger routines and/or pipelines of routines that collectively carry out different oceanographic research projects.

### SYNOPSIS

**qc_summary.R** is one component from a pipeline that collectively processes and analyzes global temperature-salinity-depth data in the ocean. This script compiles QC flags while the next script in the routine applies these flags to filter out undesired data. This was created to handle >2 million observed profiles. For more information, see the docstring in the file.

**autoarima.R** was designed to automatically and efficiently fit ARIMA models to thousands of observed ocean Lagrangian drifter velocity time series. The **R** function *auto.arima()* eliminates the need to conduct a manual parameter search for all time series and implementing it with **R**'s *apply* family of functions allows this to be done in a vectorized fashion. The parent script (not provided) compiles many such forecasts from the defined function *run.autoarima()* and writes them out to a netCDF file.

**ssa.py** demonstrates an object-oriented approach to implementing univariate and multivariate singular spectrum analysis (SSA) in Python. As noted and credited in the docstring, this started with a simple online demonstration that was motified and expanded to suit the needs of the project at hand. The working jupyter notebook shows parts of the script under development and in use.

The annotations contained in these scripts are standard practice in my coding endeavors. First, it helps me review and recall what I wrote and did when looking back after periods of time. Second, as a team player who values reproducability and sharing science, I believe any code I write for any project should be interpretable by anyone who needs to look at it. This is especially crucial for collaborative projects.

**teaching-hw-examples**: This contains examples of homework assignments I wrote and graded for an undergraduate Introduction to Physical Oceanography class in which I was tasked with teaching **R**. The students had no previous coding experience, so I incorporated it into each homework and wrote sample scripts that thoroughly explained the step-by-step processes for solving the problems. These samples were given to the students after submission. Shown here are the first and final homework assignments of the semester, demonstrating how far the students came in four months (class homework average was ~89/100). My heuristic technique with the assignments was to walk the students through a key concept in such a way that each part of the question builds upon the previous part, and then provide a practical real-world example of the concept in use. 

### PROFICIENCIES

I am proficient in both **R** and Python and generally select the language based on the task at hand: being a statistical programming language, **R** is usually best for statistical analysis and visualization and often has nifty tools that haven't been developed (or at least not released publically) in other languages, such as *auto.arima()*. But while machine learning (ML) can be done in **R** as well as in Matlab and many other languages, the ML community has largely settled on Python for its readability, simplicity, and open source status. Accordingly, I opt for Python for model development and machine learning tasks.


### CONCLUDING REMARKS

New scripts or notebooks may be added as needed or appropriate. Interests are encouraged to reach out to me with questions, requests, or solitications.

For an example of a "just for fun" side project I took up, see the respository **covid19**, which contains plots of COVID-19 data from various parts of the country (regions were selected based on the location of family and friends.) These plots are/were automatically updated daily throughout the coronovirus pandemic.


### LEGALITY NOTICE

The contents of this public respository are provided for demonstration purposes only. All rights to the products herein are reserved under standard United States copyright law. No source codes may be reproduced, distributed, or used to create derivative works without expressed permission by the owner of this respository.
