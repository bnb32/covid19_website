# README #

### What is this repository for? ###

This repo includes tools to solve a full SEAIQHRD compartmental infectious disease model. This model reduces to simpler more familiar models, such as SIR and SEIR, with the appropriate parameter values. The full model solves the following equations:

![seaiqhrd](https://latex.codecogs.com/gif.download?%5Clarge%20%5Cbegin%7Barray%7D%7Blcl%7D%20%5Cfrac%7BdS%7D%7Bdt%7D%20%26%3D%26%20-%5Cbeta%20S%20%5Cfrac%7BI+A%7D%7BN-Q-D-H%7D%20%5C%5C%20%5C%5C%20%5Cfrac%7BdE%7D%7Bdt%7D%20%26%3D%26%20%5Cbeta%20S%20%5Cfrac%7BI+A%7D%7BN-Q-D-H%7D-%5Cfrac%7BE%7D%7BT_%7Be%2Ca%7D%7D%20%5C%5C%20%5C%5C%20%5Cfrac%7BdA%7D%7Bdt%7D%20%26%3D%26%20%5Cfrac%7BE%7D%7BT_%7Be%2Ca%7D%7D-%28%5Cfrac%7Bp_%7Ba%2Ci%7D%7D%7BT_%7Ba%2Ci%7D%7D+%5Cfrac%7B1-p_%7Ba%2Ci%7D%7D%7BT_%7Ba%2Cr%7D%7D%29A%20%5C%5C%20%5C%5C%20%5Cfrac%7BdI%7D%7Bdt%7D%20%26%3D%26%20p_%7Ba%2Ci%7D%20%5Cfrac%7BA%7D%7BT_%7Ba%2Ci%7D%7D-%28%5Cfrac%7Bp_%7Bi%2Cq%7D%7D%7BT_%7Bi%2Cq%7D%7D+%5Cfrac%7Bp_%7Bi%2Cd%7D%7D%7BT_%7Bi%2Cd%7D%7D+%5Cfrac%7Bp_%7Bi%2Ch%7D%7D%7BT_%7Bi%2Ch%7D%7D+%20%5Cfrac%7B1-p_%7Bi%2Cq%7D-p_%7Bi%2Ch%7D-p_%7Bi%2Cd%7D%7D%7BT_%7Bi%2Cr%7D%7D%29I%20%5C%5C%20%5C%5C%20%5Cfrac%7BdQ%7D%7Bdt%7D%20%26%3D%26%20p_%7Bi%2Cq%7D%5Cfrac%7BI%7D%7BT_%7Bi%2Cq%7D%7D-%28%5Cfrac%7Bp_%7Bq%2Cd%7D%7D%7BT_%7Bq%2Cd%7D%7D+%5Cfrac%7Bp_%7Bq%2Ch%7D%7D%7BT_%7Bq%2Ch%7D%7D+%5Cfrac%7B1-p_%7Bq%2Ch%7D-p_%7Bq%2Cd%7D%7D%7BT_%7Bq%2Cr%7D%7D%29Q%20%5C%5C%20%5C%5C%20%5Cfrac%7BdH%7D%7Bdt%7D%20%26%3D%26%20p_%7Bq%2Ch%7D%5Cfrac%7BQ%7D%7BT_%7Bq%2Ch%7D%7D+p_%7Bi%2Ch%7D%5Cfrac%7BI%7D%7BT_%7Bi%2Ch%7D%7D-%28%5Cfrac%7Bp_%7Bh%2Cd%7D%7D%7BT_%7Bh%2Cd%7D%7D+%5Cfrac%7B1-p_%7Bh%2Cd%7D%7D%7BT_%7Bh%2Cr%7D%7D%29H%20%5C%5C%20%5C%5C%20%5Cfrac%7BdR%7D%7Bdt%7D%20%26%3D%26%20%281-p_%7Bq%2Ch%7D-p_%7Bq%2Cd%7D%29%5Cfrac%7BQ%7D%7BT_%7Bq%2Cr%7D%7D%20+%20%281-p_%7Bi%2Cq%7D-p_%7Bi%2Ch%7D-p_%7Bi%2Cd%7D%29%5Cfrac%7BI%7D%7BT_%7Bi%2Cr%7D%7D%20+%20%281-p_%7Ba%2Ci%7D%29%5Cfrac%7BA%7D%7BT_%7Ba%2Cr%7D%7D+p_%7Bh%2Cr%7D%5Cfrac%7BH%7D%7BT_%7Bh%2Cr%7D%7D%20%5C%5C%20%5C%5C%20%5Cfrac%7BdD%7D%7Bdt%7D%20%26%3D%26%20p_%7Bq%2Cd%7D%20%5Cfrac%7BQ%7D%7BT_%7Bq%2Cd%7D%7D%20+%20p_%7Bi%2Cd%7D%5Cfrac%7BI%7D%7BT_%7Bi%2Cd%7D%7D+p_%7Bh%2Cd%7D%5Cfrac%7BH%7D%7BT_%7Bh%2Cd%7D%7D%20%5C%5C%20%5Cend%7Barray%7D)

The S,E,A,I,Q,H,R,D compartments represent Susceptible, Exposed, Asymptomatic-Infected, Symptomatic-Infected, Self-Quarantined, Hospitalized, Recovered, and Dead members of the total population N. 

All parameters of the following form represent mean times for transitions between compartments:

![params](https://latex.codecogs.com/gif.download?%5Clarge%20T_%7Bx%2Cy%7D%20%3A%20X%20%5Crightarrow%20Y)

Parameters of the following form represent probabilities for transitions between compartments:

![probs](https://latex.codecogs.com/gif.download?%5Clarge%20p_%7Bx%2Cy%7D%20%3A%20X%20%5Crightarrow%20Y)


### How do I get set up? ###

To solve the previous equations initial values for all compartments and parameter values must be specified. The repo is currently setup to download `new_cases`, `total_cases`, `deaths`, and `hospitalizations` from COVID case sources, in addition to population data `(N)` for a specified state. Using this data compartments can be initialized and probabilities for death and hositalization can be computed. Transition times between compartments must be specified manually, as well as the probabilities for some transitions. 

**Current COVID research suggests a recovery time of 14 days `(Tir: I->R)`, an incubation time of 5 days `(Tai: A->I)`, a time from development of symptoms to death 30 days `(Tid: I->D)`, and a probability of developing symptoms once infected of roughly 0.6 `(pai: A->I)`. As COVID does not seem to have a non-contagious latency period we use a minimum time from exposure to asymptomatic `(Tea: E->A)` of one time step. We also assume once an individual recognizes symptoms they will self-quarantine in a single time step `(Tiq: I->Q)` with a probability of 0.9 `(piq: I->Q)`.**  

The doubling rate `(Td)`, and subsequently `beta`, can be computed by approximating initial growth as exponential and performing a linear regression:

![doubling](https://latex.codecogs.com/gif.download?%5Clarge%20%5Cbegin%7Barray%7D%7Blcl%7D%20I_%7Btotal%7D%20%26%20%5Capprox%20%26%20I_%7B0%7D2%5E%7Bt/T_%7Bd%7D%7D%5Capprox%20I_%7B0%7D2%5E%7Bt%28%5Cbeta%20N/S_%7B0%7D%29%7D%20%5C%5C%20%5CDelta%20I_%7Btotal%7D%20%26%20%5Capprox%20%26%20%5Cfrac%7Bln%282%29%7D%7BT_%7Bd%7D%7DI_%7Btotal%7D%20%5CDelta%20t%20%5C%5C%20%5Cfrac%7Bln%282%29%7D%7BT_%7Bd%7D%7D%20%26%20%5Capprox%20%26%20%5Cbeta%20%5Cfrac%7BN%7D%7BS_%7B0%7D%7D%20%5Cend%7Barray%7D)

The length of data to use for this regression can be manually specified. 

The main routine in `__init__.py` takes the following arguments:

`state`: US state  
`county`: US county  
`n_days`: number of days to simulate  
`Tea`: mean time from exposure to asymptomatic (latency time)  
`Tiq`: mean time from recognizing symptoms to self-quarantine  
`Tai`: mean time to develop symptoms from initial infection (incubation time)  
`Tir`: mean time to recover after developing symptoms  
`Tid`: mean time from developing symptoms to death  
`Tih`: mean time to be hospitalized after developing symptoms  
`Tqh`: mean time to be hospitalized after self-quarantine  
`piq`: probability to self-quarantine after recognizing symptoms  
`pai`: probability to develop symptoms after infection (symptomaticity ratio)  
`Sd`: fractional reduction in physical contact  
`Sd_period`: duration of contact reduction  
`Sd_delay`: time to implement contact reduction  
`detection_rate`: probability of detecting an infection  
`data_days`: number of days of data to use for computing doubling time `(Td)`  
`n_substeps`: time steps per day  
`R0`: reproductive rate (calculate from data if negative)  

* All models are solved using 4th order Runge-Kutta, the details of which are in `models.py`. Functions for downloading data and initializing model are in `fetch.py`. Functions for plotting are in `display.py`. The main routine to run the model and plot is in `__init__.py`.
* Required packages include: matplotlib, numpy, scipy

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
