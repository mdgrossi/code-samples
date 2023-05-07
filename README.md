# code-samples

#####Sample programming scripts and other technical capabilities to demonstrate data processing, analysis, and management proficiencies

## SUMMARY

These scripts are intended to be short samples that demonstrate my capabilities as a scientist at handling, processing, and analyzing large data sets using various statistical methods pertinent to different research projects I have conducted over the years. Each script represents a glimpse into much larger routines and/or pipelines of routines that collectively carry out different oceanographic research projects.

This repository, whose content is described below, is organized as follows:

```bash
.
├── python_examples
│   ├── ssa.py
│   ├── singularSpectrumAnalysis.ipynb
├── r_examples
│   ├── autoarima.R
│   ├── qc_summary.R
├── teaching_hw_examples
│   ├── EOS80.R
│   ├── HW01_sample_code.R
│   ├── HW01_GoM_data.csv
│   ├── HW01_answers.pdf
│   ├── HW05_sample_code.R
│   ├── HW05_answers.pdf
├── README.md
└── .gitignore
```

## PROFICIENCIES

I am proficient in both **R** and Python and generally select the language based on the task at hand: being a statistical programming language, **R** is usually best for statistical analysis and visualization and often has nifty tools that haven't been developed (or at least not released publically) in other languages, such as *auto.arima()*. But while machine learning (ML) can be done in **R** as well as in Matlab and many other languages, the ML community has largely settled on Python for its readability, simplicity, and open source status. Accordingly, I opt for Python for model development and machine learning tasks. Example scripts from both languages are provided.

## /docker\_examples

_Coming soon_

## /python\_examples

<details>
    <summary>
        Python scripts and Jupyter notebooks
    </summary>
    
**ssa.py** in demonstrates an object-oriented approach to implementing univariate and multivariate singular spectrum analysis (SSA) in Python. As noted and credited in the docstring, this started with a simple online demonstration that was modified and expanded to suit the needs of the project at hand.

**singularSpectrumAnalysis.ipynb** is a scratch Jupyter notebook showing parts of the **ssa.py** script under development and in use. It is used for testing, experimentation, plotting, etc. This is a representative "working notebook"; it is not intended to a something one might use to showcase, narrate, or inform about the development of some process.

</details>

### /r\_examples

<details>
    <summary>
        R scripts
    </summary>

**qc_summary.R** is one component from a pipeline that collectively processes and analyzes global temperature-salinity-depth data in the ocean. This script compiles QC flags while the next script in the routine applies these flags to filter out undesired data. This was created to handle >2 million observed profiles. For more information, see the docstring in the file.

**autoarima.R** was designed to automatically and efficiently fit ARIMA models to thousands of observed ocean Lagrangian drifter velocity time series. The **R** function *auto.arima()* eliminates the need to conduct a manual parameter search for all time series and implementing it with **R**'s *apply* family of functions allows this to be done in a vectorized fashion. The parent script (not provided) compiles many such forecasts from the defined function *run.autoarima()* and writes them out to a netCDF file.

</details>

## /teaching\_hw\_examples

<details>
    <summary>
        Teaching heuristics and homework assignments: sample question sets with answers, scripts to solve the problems, and data files where appropriate
    </summary>

**HW01\_\*** and **HW05\_\***: Examples of homework assignments I wrote and graded for an undergraduate Introduction to Physical Oceanography class in which I was tasked with teaching **R**. The students had no previous coding experience, so I incorporated it into each homework and wrote sample scripts that thoroughly explained the step-by-step processes for solving the problems. These samples were given to the students after the homework assignment was due. Shown here are the first and final homework assignments of the semester, demonstrating how far the students came in four months (class homework average was ~89/100). My heuristic technique with the assignments was to walk the students through a key concept in such a way that each part of the question builds upon the previous part, and then provide a practical real-world example of the concept in use. The question sets with answers are provided for context as well as the data used in HW01.

**EOS80.R**: Simple Equation of State (EOS) of Seawater (1980) function I created and provided to the students for use in HW01.

</details>


## /latex\_examples

_Coming soon_

### CONCLUDING REMARKS

The annotations contained in these scripts are standard practice in my coding endeavors. First, it helps me review and recall what I wrote and did when looking back after periods of time. Second, as a team player who values reproducability and sharing science, I believe any code I write for any project should be interpretable by anyone who needs to look at it. This is especially crucial for collaborative projects.

New scripts or notebooks may be added to this repository as needed. Interests are encouraged to reach out to me with questions, requests, or solitications.

For an example of a "just for fun" side project I took up, see the respository **covid19**, which contains plots of COVID-19 data from various parts of the country (regions were selected based on the location of family and friends.) These plots were automatically updated daily throughout the pandemic. Daily updating was terminating when the data sources themselves ceased providing reliable or timely data.


### LEGALITY NOTICE

_The contents of this public respository are provided for demonstration purposes only. All rights to the products herein are reserved under standard United States copyright law. No source codes may be reproduced, distributed, or used to create derivative works without expressed permission by the owner of this respository._
