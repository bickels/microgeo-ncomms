# microGEO-ncomms

μGEO - micro geography of microbial abundance and diversity across spatial scales. Minimal example demonstrating the model implementation.
---

## System requirements
### Tested on: 
Windows 7 64bit
### Dependencies (tested version):
	- python (2.7.14)
		- numpy (1.14.3)
		- scipy (1.1.0)
		- pandas (0.22.0)
		- matplotlib (2.2.2)
## Installation
Once python and dependencies are installed the script (`example.py`) can be executed from the command line after navigating to the same directory:
```
C:\> cd path\to\script\
C:\path\to\script\python example.py
```
Installation (Anaconda distribution recommended) should be possible within less than 30min.

## Demo
The demonstration includes a small excerpt of data from the Earth Microbiome Project as published by Thompson *et al.* 2017 (DOI:10.1038/nature24621) combined with covariates as described in SI Table S1 of the submitted manuscript. 
Running the demo `example.py` produces a figure comparing two modeled species abundance distributions (one and multiple species per habitat) with the observed distribution using the georeferenced metadata of the sample. It additionally prints the samples bacterial species richness and evenness. The runtime should be less then 5s.
