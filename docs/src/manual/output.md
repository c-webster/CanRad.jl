

## Outputs

### Sky-view fraction (svf)

**svf plane:** Sky-view fraction calcuated for the perspective of a horizontal flat uplooking surface (sensor or ground). Calculated by weighting zenith rings in synthetic hemispheric image according to their surface area projected onto a flat horizontal surface

**svf hemi:** Sky-view fraction calculated for the perspective of a hemispherically shaped surface (sensor or e.g. plant). Calculated by weighting zenith rings in synthetic hemispheric image according to their surface area on the hemisphere

### Forest transmissivity

Time-varying direct beam shortwave transmissivity. Calculated by determining ratio of canopy/sky pixels in front of the projected position of the solar disc.

Calculated at 2 minute intervals to ensure a temporally complete solar track, then averaged to the user-defined interval (*tstep*). The output time stamps are the beginning of the averged period. 



### SWR radiation

**SWR_total**
Total incoming shortwave radiation (diffuse + direct). 

**SWR_direct**
Direct incoming shortwave radiation in the path of the solar position.

**SWR_diffuse**
Diffuse incoming shortwave radiation calculated from sky-view fraction


*Calculated within T/L/C2R:*
- calculated at the same 2 minute interval as forest transmissivity, then averaged
- SWR_diffuse is not saved as output; can be calculated using the difference between total and direct radiation. 
- only calculates maximum potential value, using an atmospheric transmissivity of 1. 

*Calculated with CalcSWRFromNC.jl*
- calculated at the temporal resolution of the saved forest transmissivity variable.
- can calculate real incoming shortwave radiation if above-canopy data is available. !! note this data must match the time steps of the transmissivity data in the .nc file exactly.



### Synthetic hemispheric images

The settings file has an option to save the calculated synthetic images in a netcdf file. These are created as .png files in a separate folder. 
