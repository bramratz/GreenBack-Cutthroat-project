# Estimating Ancestry of F2 Hybrid Offspring of Greenback Cutthroat Trout and Colorado River Trout Parents

The Greenback Cutthroat trout, a subspecies of Colorado River Cutthroat trout, was believed to be extinct in its native range by 1937. Its subsequent rediscovery in the 1970s have spurned an ongoing convervation program in Colorado aimed at propogating and increasing genetically pure populations of the species. However, previous attempts to propogate the species in hacheries have been unsuccessful, presumably due to high genetic loads and inbreeding depression. As a result individuals propogated in hatcheries will likely never be released back into their native range. In collaboration with Colorado Parks and Wildlife 24 hybrid families of Greenback Cutthroat and Colorado River Cutthroat have been propagated to the F2 generation. Using these families we intend to estimate the overall ancestry in each parental subspecies for individual offspring and quantify variation across species. It is expected that introducing Colorado River Cutthroat to an inbred population of Greenback Cutthroat will artificially create gene flow into the population, resulting in hybrid offspring with higher fitness relative to pure offspring.

## Goals
1. Estimate the overall ancestry in each parental subspecies for individual offspring and quantify variation within and across families 
2. Quantify locus-specific ancestry in hybrids to identify which alleles are over/under-represented in hybrids relativeto expectations. Compare against Rainbow trout (RBT) genome annotations for hints at potential function

## Main computational tasks 
1. Bash filtering script
  - Filters vcfs created using **smatools** for missing data at the locus and individual level, minor allele frequency, and paralogs
  - Creates output file in .mpgl format to used by entropy
2. estimate ancestry of individuals using **entropy**

