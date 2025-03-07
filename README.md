# X-ray albedo corrections for solar flare Hard X-ray Spectroscopy

X-ray albedo green matrix correction for solar hard X-ray analysis 

Ray tracing and scattering based on the following paper:
https://ui.adsabs.harvard.edu/abs/2006A%26A...446.1157K/abstract

To correct for albedo in OSPEX see the manual:
https://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_explanation.htm#Albedo%20Correction

X-ray albedo included in GX simulator: 
https://github.com/Gelu-Nita/GX_SIMULATOR/blob/d772a411705ceed9db2a5b6eadfa14cd5dda190f/userslib/xray/xray_tt_albedo.pro


To generate Green Matrices, from the IDL command prompt (SSW IDL), run the following:

    IDL> make_all_green_1kev

The routine will generate IDL save files in the local folder:
> ls                                             
green_compton_mu005.dat  green_compton_mu055.dat 
green_compton_mu010.dat  green_compton_mu060.dat 
green_compton_mu015.dat  green_compton_mu065.dat 
green_compton_mu020.dat  green_compton_mu070.dat 
green_compton_mu025.dat  green_compton_mu075.dat 
green_compton_mu030.dat  green_compton_mu080.dat 
green_compton_mu035.dat  green_compton_mu085.dat 
green_compton_mu040.dat  green_compton_mu090.dat 
green_compton_mu045.dat  green_compton_mu095.dat 
green_compton_mu050.dat                          

and PostScript reference albedo spectrum:
> ls     
albedo_mu077_pow2_reference.ps

See discussion here: https://ui.adsabs.harvard.edu/abs/2007A%26A...466..705K/abstract



