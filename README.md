# klein2013-cloud-error-metrics
Compute four scalar measures of the fidelity of model cloud simulations described in Section 4 of Klein et al. (2013). 

References
----------
Klein, S.A., Y. Zhang, M.D. Zelinka, R.N. Pincus, J.Boyle, and P.J. Gleckler, 2013: Are climate model simulations of clouds improving? An evaluation using the ISCCP simulator. J. Geophys. Res. 118, 1329-1342. doi: 10.1002/jgrd.50141.

Zelinka, M. D., S. A. Klein, and D. L. Hartmann, 2012: Computing and Partitioning Cloud Feedbacks Using Cloud Property Histograms. Part I: Cloud Radiative Kernels. J. Climate, 25, 3715–3735. doi:10.1175/JCLI-D-11-00248.1.

Input
----------

The code makes use of the following data, all of which are available at https://github.com/mzelinka/klein2013-cloud-error-metrics/tree/master/data:


| Frequency | Name | Description | Unit | File Format |
|:----------|:-----------------------------|:-------------|:------|:------------|
| monthly mean | AC_clisccp | climatological mean annual cycle of observed ISCCP cloud fraction histograms | % | nc | 
| monthly mean | AC_clmodis | climatological mean annual cycle of observed MODIS cloud fraction histograms | % | nc |
| monthly mean | clisccp | model ISCCP simulator cloud fraction histograms | % | nc |
| monthly mean | rsuscs | model upwelling SW flux at the surface under clear skies | W/m^2 | nc |
| monthly mean | rsdscs | model downwelling SW flux at the surface under clear skies | W/m^2 | nc |
| monthly mean | LWkernel | LW cloud radiative kernel | W/m^2/% | nc |
| monthly mean | SWkernel | SW cloud radiative kernel | W/m^2/% | nc |


The GCM simulator-oriented ISCCP and MODIS cloud observational benchmarks from which AC_clisccp and AC_clmodis are derived are available on the CFMIP-OBS page http://climserv.ipsl.polytechnique.fr/cfmip-obs/

Output
----------
The following 4 scalars:

E_TCA: the total cloud amount error 

E_CP: the cloud properties error 

E_LWCP: the LW-relevant cloud properties error 

E_SWCP: the SW-relevant cloud properties error 


Figure Generation
----------
Is a script to draw a figure in the paper included? No
