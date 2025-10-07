# Clinical Trial Side-Effects Analysis: A Multilevel Modeling Approach
Clinical Trial Side-Effects Analysis: A Multilevel Modeling Approach

[![View Live Report](https://img.shields.io/badge/View-Live%20Report-brightgreen)](https://htmlpreview.github.io/?https://github.com/Vismaya-SD/Trial-Side-Effects-Analysis-Three-Level-Linear-Mixed-Effects-Model/blob/main/Trial-Side-Effects-Analysis.html)
[![Language](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)

## Project Overview

This repository contains a comprehensive analysis of simulated data from a longitudinal clinical trial. The project's goal is to model the severity of drug-induced side effects over time while accounting for the complex, hierarchical structure of the data (i.e., multiple measurements nested within patients, who are in turn nested within different clinical sites).

The analysis demonstrates the application of mixed-effects models, a powerful statistical technique for handling such correlated data structures common in medical research.

## Key Methods Used

* **Three-Level Linear Mixed-Effects Model**: Used to analyze the continuous `side_effect_score` (0-10 scale). This was the primary and most informative model.
* **Generalized Linear Mixed-Effects Model**: Used for a secondary analysis on a binary outcome (`side_effect_present`), demonstrating how to handle logistic regression within a mixed-effects framework.
* **Exploratory Data Analysis (EDA)**: Data visualization to understand trends and distributions before modeling.
* **Likelihood Ratio Tests (LRT)**: To formally test the significance of the random effects.

## Repository Structure

* `Trial-Side-Effects-Analysis.Rmd`: The R Markdown source file. It contains all the R code for data simulation, EDA, statistical modeling, and the complete narrative of the analysis.
* `Trial-Side-Effects-Analysis.html`: The final, self-contained HTML report generated from the R Markdown file. **Click the "View Live Report" badge at the top to see the rendered page directly.**
* `.gitignore`: A standard file that tells Git to ignore temporary files specific to R and RStudio.

## How to Reproduce the Analysis

To reproduce this analysis on your own machine, follow these steps:

1.  Clone or download this repository.
2.  Open the `Trial-Side-Effects-Analysis.Rmd` file in RStudio.
3.  Ensure you have the required R packages installed by running the following command in your R console:
    ```r
    install.packages(c("tidyverse", "lme4", "lmerTest", "sjPlot", "patchwork", "knitr"))
    ```
4.  Click the **Knit** button in RStudio to re-generate the `Trial-Side-Effects-Analysis.html` report.

## Summary of Findings

The analysis revealed several key insights:

1.  **Continuous Score is More Informative**: The linear mixed-effects model on the continuous `side_effect_score` was highly effective, identifying significant effects for `time`, `drug`, `age`, and `sex`.
2.  **Information Loss in Binary Model**: Converting the outcome to a binary variable (`present`/`absent`) resulted in a significant loss of statistical power. The binary model was unstable (failed to converge properly) and only found `time` to be a significant predictor.
3.  **Importance of Hierarchical Modeling**: Both `site` and `patient` levels contributed significantly to the variance in the data (accounting for over 27% of total variance combined), confirming that a mixed-effects model was the appropriate choice over a standard linear regression.
