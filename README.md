# VVO-project-MGarcia - Final Degree Project Scripts
Here I present the scripts developed throughout my Final Degree Project (TFG). 

My project has focused on the optimization and robustness testing of **VHIO’s Visual Omics (VVO)**, an interactive application under development. Based on `R Shiny` and `iSEE`, VVO is designed for the visual exploration of bulk transcriptomic data obtained through the internal workflows of the **VHIO Bioinformatics Unit**.

---

## Project Objectives

* **Feature Optimization:** I took over the project from a previous version and my main aim was to optimize it by adding two new fields of analysis: Differential Expression Analysis (**DEA**) and Functional Enrichment Analysis (**FEA**).
* **Robustness Testing:** I tested the application's robustness by uploading different datasets transformed to accommodate the input format required by VVO (a `DeeDeeExperiment` object).

---

## File Descriptions

The repository contains the following developed files:

* **`VVO.R`**: The main script containing the interactive application code and the newly implemented features.
* **`sceCreation.R`** & **`ddeCreations.R`**: Scripts specifically developed to transform and format different datasets into `DeeDeeExperiment` objects for robustness testing.
