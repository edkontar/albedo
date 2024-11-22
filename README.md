# albedo

X-ray albedo green matrix correction for solar hard X-ray analysis 

Ray tracing and scattering based on the following paper:
https://ui.adsabs.harvard.edu/abs/2006A%26A...446.1157K/abstract

To correct for albedo in OSPEX see the manual:
https://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_explanation.htm#Albedo%20Correction

X-ray albedo included in GX simulator: 
https://github.com/Gelu-Nita/GX_SIMULATOR/blob/d772a411705ceed9db2a5b6eadfa14cd5dda190f/userslib/xray/xray_tt_albedo.pro


To generate Green Matrices, from the IDL command prompt (SSW IDL), run the following:

    IDL> make_all_green_1kev

The routine will generate matrices in the local folder.



