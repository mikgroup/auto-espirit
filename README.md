# Auto-ESPIRiT demonstration code.

Code used to generate figures for the MRM manuscript,
__SURE-based Automatic Parameter Selection For ESPIRiT Calibration__.

Written by Siddharth Iyer. Please feel free to post an issue on the repository page if there is a problem.  

This code may be freely used and modified for educational, research, and not-for-profit purposes (See
[`LICENSE`](LICENSE) for more information).

## Structure

- `dem_code` contains code used to generate the data in Figure 1.
- `exp_code` contains code used to generate the data in Figure 4.
- `fig_code` consists for helper functions to generate the figures themselves.

## Acknowledgements

Divergence code in `utils/svWeightsSure.m` was obtained and modified from
http://statweb.stanford.edu/~candes/SURE/data.html. It was originally
written by E.J. Candès, C.A. Sing-Long, and J.D. Trzasko.

Some Matlab utility functions were written by Michael Lustig in the ESPIRiT
Matlab reference implementation. They can be obtained from
http://people.eecs.berkeley.edu/~mlustig/Software.html.

All rights/distribution are the same as for the original code,
and should cite the original author and webpage.

## DOI
[![DOI](https://zenodo.org/badge/145603994.svg)](https://zenodo.org/badge/latestdoi/145603994)
