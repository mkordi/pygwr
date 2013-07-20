Tokyo Mortality example for pygwr
=================================

This example shows how to use pygwr software.

The example uses the Tokyo mortality dataset of Nakaya et al. 2005, and is available from this website: http://www.st-andrews.ac.uk/geoinformatics/gwr/gwr-downloads/
Only the file Tokyomortality.txt is needed.

Once the file is downloaded and placed in the example folder, it is possible to run the GWR by running the tokyo_gwr.py script.

The example runs a model similar to the one provided with GWR4, but with some simplification, because pygwr does not provide the possibility to include global variables. The provided example calculates a Poisson GWR with db2564 as dependent variable, and OCC_TEC and UNEMP as independent variables. A fixed Gaussian kernel is used, and for bandwidth selection, a simple interval search between 5000 and 20000 with steps of 1000 is performed.

