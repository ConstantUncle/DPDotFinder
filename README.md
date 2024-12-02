Simple data analysis tool made in Python as a part of internship with JAXA. Takes X-Ray telescope event file, estimated period, and estimated first time derivative of the period, and iterates over the EFSearch tool (part of [FTOOLS](https://heasarc.gsfc.nasa.gov/ftools/), a suite of data analysis tools for FITS files) multiple times to find best estimate for both. Output is in the form of a contour heat map, showing chi squared values plotted against Period and PDot.

Requires FTOOLS
