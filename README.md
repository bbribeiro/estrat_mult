# Endogenous stratification in RCTs with multiple treatments

Author: Bernardo Ribeiro

Last update: December 1st, 2018

This code is an extension of Stata package `estrat`, written by Jeremy Ferwerda. The original package `estrat` implements the algrothim proposed by Abadie, Chingos and West (2018) to calculate heterogeneous treatment effects in RCTs based on baseline characteristics that predict outcomes of interest. The marginal changes in `estrat_mult` allow for

- multiple treatment arms

- estratification based on (potential) outcomes that are not necessarily the same as the final outcomes of interest

### Syntax

`estrat depvar treatment [if] [in] [, options]`

Since the variable used in the estratification may differ from the outcome of interest, we now have to specify both in the options.

  - `depvar` specifies the final outcome of interest, as in `estrat`

  - `treatment` specifies the list of variable with the treatment assignment, with multiple treatments allowed. In this case, `treatment' should be a list of binary indicators, one for each treatment arm not including control arm

The command requires the following options

  - `potential` specifies the (potential) outcome to be predicted and used in the estratification

  - `pred` specifies the variables used to generate predicted outcomes (i.e. `potential`)

  - `control` specifies the control group (binary indicator)
  
Following the original program, the following are optional

  - `groups` specifies the number of sub-groups
  
  - `cov` specifies covariates used when estimating treatment effects
  
  - `reps` specifies the number of repeated split sample (RSS) repetitions
  
  - `boot` specifies the number of bootstrap repetitions
  
  - `loo_only` instructs the program to only estimate leave-one-out (LOO) results
  
  - `rss_only` instructs the program to only estimate repeated-split-sample (RSS) results
  
  - `savegroup` saves the group number assigned to each observation in the LOO calculation
  
### Codes

The folder 'ado' contains the ado-file of the original package `estrat` and the modified version `estrat_mult`

### Example

The file 'jtpa.do' is a do-file that runs the example in the original help file of `estrat` and modifies the data to include multiple treatments and a random variable to be used as potential outcome (which is not the final outcome of interest) and runs `estrat_mult`

### References

Abadie, Alberto, Matthew Chingos and Martin West, 2018. "Endogenous Stratification in Randomized Experiments," The Review of Economics and Statistics, vol 100(4), pages 567-580
