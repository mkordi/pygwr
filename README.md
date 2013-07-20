pygwr: a simple GWR in Python
=============================

Python library for Geographically Weighted Regression. Both Gaussian and Poisson GWR are supported.

pygwr builds on top of the statsmodels Python package (http://statsmodels.sourceforge.net). statsmodels provides all statistical algorithms underlying to GWR. pygwr uses a slightly modified version of statsmodels for supporting geographically weighted Poisson regression. pygwr implements all the weighting scheme of GWR.

An example of how to use pygwr can be found in the example folder, in the tokyo_gwr.py script. The readme.txt file in the example folder gives the instructions where to find the required dataset.


Copyright
---------

Copyright 2013 Maryam Kordi

pygwr is distributed under an open-source licence, more specifically the GNU Public Licence (GPL) version 3 or later.

The modified statsmodels package within pygwr (gwstatsmodels) is under its initial licence, which is a modified BSD licence. See the statsmodels website at http://statsmodels.sourceforge.net for more information.


Licence (GPL)
-------------

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.


