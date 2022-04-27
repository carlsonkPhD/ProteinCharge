# ProteinCharge

R package to calculate the theoretical charge of a protein at different pH values. 

## Installation

To install ProteinCharge directly from github, the devtools package is required. Run the following to install devtools:

```R
install.packages("devtools")
```

Once devtools is installed, ProteinCharge can be installed with the following:

```R
library("devtools")
install_github("calrsonkPhD/ProteinCharge", dependencies = TRUE)
```

## Loading

To load ProteinCharge use the following:

```R
library(ProteinCharge)
```

## Usage

The ProteinCharge package contains a single function named 'titration'. The help entry for this function can be accessed with

```R
?titration
```

The amino acid sequence of the protein can be entered either as a string passed to the 'sequence' argument or in the form of a txt file passed to the 'read_source' argument. The theoretical charge of the protein at each point along the titration is printed to StdOut and the graph of these values appears in the Plots pane of RStudio.