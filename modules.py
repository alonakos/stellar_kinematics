import astropy.units as u
import numpy
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
import pandas as pd

def get_stars (name, ra, dec, radius = 0.0194):
    
        #create a query for starts with ra, dec 
        #,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val 
#         query = f"SELECT TOP 500 gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.pmra,gaia_source.pmdec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.dr2_radial_velocity,gaia_source.dr2_radial_velocity_error \
#             FROM gaiaedr3.gaia_source \
#             WHERE \
#             gaia_source.dr2_radial_velocity IS NOT NULL \
#             AND \
#             gaia_source.parallax IS NOT NULL \
#             AND \
#             CONTAINS( \
#             POINT('ICRS', gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
#             CIRCLE('ICRS', {ra}, {dec} ,{radius}) \
#             )=1" 
        query = f"SELECT TOP 500 gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.pmra_pmdec_corr,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.dr2_radial_velocity,gaia_source.dr2_radial_velocity_error \
             FROM gaiaedr3.gaia_source \
             WHERE \
             gaia_source.dr2_radial_velocity IS NOT NULL \
             AND \
             gaia_source.parallax IS NOT NULL \
             AND \
             CONTAINS( \
             POINT('ICRS', gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
             CIRCLE('ICRS', {ra}, {dec} ,{radius}) \
             )=1"
 
        # query = f"SELECT TOP 500 gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.pmra_pmdec_corr,gaia_source.ruwe,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.non_single_star,gaia_source.has_xp_continuous,gaia_source.has_xp_sampled,gaia_source.has_rvs,gaia_source.has_epoch_photometry,gaia_source.has_epoch_rv,gaia_source.has_mcmc_gspphot,gaia_source.has_mcmc_msc,gaia_source.teff_gspphot,gaia_source.logg_gspphot,gaia_source.mh_gspphot,gaia_source.distance_gspphot,gaia_source.distance_gspphot_lower,gaia_source.distance_gspphot_upper,gaia_source.azero_gspphot,gaia_source.ag_gspphot,gaia_source.ebpminrp_gspphot, \
        #     FROM gaiaedr3.gaia_source \
        #     WHERE \
        #     gaiaedr3.gaia_source.dr2_radial_velocity IS NOT NULL \
        #     AND \
        #     gaia_source.parallax IS NOT NULL \
        #     AND \
        #     CONTAINS( \
        #     POINT('ICRS', gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
        #     CIRCLE('ICRS', {ra}, {dec} ,{radius}) \
        #     )=1"
        
        #launch a query

        job = Gaia.launch_job_async (query, output_file="star_data/neighbouring_stars/"+name+".vot", dump_to_file=True)

        return job.get_results()

def get_halo_stars (p_e):
    
        #create a query for starts with parallax/parallax error is more than 5
        query = f"SELECT gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.pmra,gaia_source.pmdec, gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.dr2_radial_velocity,gaia_source.dr2_radial_velocity_error \
            FROM gaiaedr3.gaia_source  \
            WHERE \
             gaiaedr3.gaia_source.parallax / gaiaedr3.gaia_source.parallax_error > {p_e} \
             AND \
             gaia_source.dr2_radial_velocity IS NOT NULL" 
        
        #launch a query
        
        job = Gaia.launch_job_async (query, dump_to_file=True)

        return job.get_results()
    

def gal_uvw(distance=None, lsr=None, ra=None, dec=None, pmra=None, pmdec=None, vrad=None, plx=None):
   """
    NAME:
        GAL_UVW
    PURPOSE:
        Calculate the Galactic space velocity (U,V,W) of star
    EXPLANATION:
        Calculates the Galactic space velocity U, V, W of star given its
        (1) coordinates, (2) proper motion, (3) distance (or parallax), and
        (4) radial velocity.
    CALLING SEQUENCE:
        GAL_UVW [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD= , DISTANCE=
                 PLX= ]
    OUTPUT PARAMETERS:
         U - Velocity (km/s) positive toward the Galactic *anti*center
         V - Velocity (km/s) positive in the direction of Galactic rotation
         W - Velocity (km/s) positive toward the North Galactic Pole
    REQUIRED INPUT KEYWORDS:
         User must supply a position, proper motion,radial velocity and distance
         (or parallax).    Either scalars or vectors can be supplied.
        (1) Position:
         RA - Right Ascension in *Degrees*
         Dec - Declination in *Degrees*
        (2) Proper Motion
         PMRA = Proper motion in RA in arc units (typically milli-arcseconds/yr)
         PMDEC = Proper motion in Declination (typically mas/yr)
        (3) Radial Velocity
         VRAD = radial velocity in km/s
        (4) Distance or Parallax
         DISTANCE - distance in parsecs
                    or
         PLX - parallax with same distance units as proper motion measurements
               typically milliarcseconds (mas)
   
    OPTIONAL INPUT KEYWORD:
         /LSR - If this keyword is set, then the output velocities will be
                corrected for the solar motion (U,V,W)_Sun = (-8.5, 13.38, 6.49)
                (Coskunoglu et al. 2011 MNRAS) to the local standard of rest.
                Note that the value of the solar motion through the LSR remains
                poorly determined.
     EXAMPLE:
         (1) Compute the U,V,W coordinates for the halo star HD 6755.
             Use values from Hipparcos catalog, and correct to the LSR
         ra = ten(1,9,42.3)*15.    & dec = ten(61,32,49.5)
         pmra = 627.89  &  pmdec = 77.84         ;mas/yr
         dis = 144    &  vrad = -321.4
         gal_uvw,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,dis=dis,/lsr
             ===>  u=154  v = -493  w = 97        ;km/s
   
         (2) Use the Hipparcos Input and Output Catalog IDL databases (see
         http://idlastro.gsfc.nasa.gov/ftp/zdbase/) to obtain space velocities
         for all stars within 10 pc with radial velocities > 10 km/s
   
         dbopen,'hipparcos,hic'      ;Need Hipparcos output and input catalogs
         list = dbfind('plx>100,vrad>10')      ;Plx > 100 mas, Vrad > 10 km/s
         dbext,list,'pmra,pmdec,vrad,ra,dec,plx',pmra,pmdec,vrad,ra,dec,plx
         ra = ra*15.                 ;Need right ascension in degrees
         GAL_UVW,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,plx = plx
         forprint,u,v,w              ;Display results
    METHOD:
         Follows the general outline of Johnson & Soderblom (1987, AJ, 93,864)
         except that U is positive outward toward the Galactic *anti*center, and
         the J2000 transformation matrix to Galactic coordinates is taken from
         the introduction to the Hipparcos catalog.
    REVISION HISTORY:
         Written, W. Landsman                       December   2000
         fix the bug occuring if the input arrays are longer than 32767
           and update the Sun velocity           Sergey Koposov June 2008
   	   vectorization of the loop -- performance on large arrays
           is now 10 times higher                Sergey Koposov December 2008
   """
   
   n_params = 3
   
   if n_params == 0:   
      print( 'Syntax - GAL_UVW, U, V, W, [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD=')
      print( '                  Distance=, PLX=')
      print( '         U, V, W - output Galactic space velocities (km/s)')
      return None
   
   if ra is None or dec is None:   
      raise Exception('ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)')
   if plx is None and distance is None:
      raise Exception('ERROR - Either a parallax or distance must be specified')
   if distance is not None:
      if numpy.any(distance==0):
         raise Exception('ERROR - All distances must be > 0')
      plx = 1e3 / distance          #Parallax in milli-arcseconds
   if plx is not None and numpy.any(plx==0):   
      raise Exception('ERROR - Parallaxes must be > 0')
   
   cosd = numpy.cos(numpy.deg2rad(dec))
   sind = numpy.sin(numpy.deg2rad(dec))
   cosa = numpy.cos(numpy.deg2rad(ra))
   sina = numpy.sin(numpy.deg2rad(ra))
   
   k = 4.74047     #Equivalent of 1 A.U/yr in km/s   
   a_g = numpy.array([[0.0548755604, +0.4941094279, -0.8676661490],
                [0.8734370902, -0.4448296300, -0.1980763734], 
                [0.4838350155, 0.7469822445, +0.4559837762]])
   
   vec1 = vrad
   vec2 = k * pmra / plx
   vec3 = k * pmdec / plx
   
   u = (a_g[0,0] * cosa * cosd + a_g[1,0] * sina * cosd + a_g[2,0] * sind) * vec1 + (-a_g[0,0] * sina + a_g[1,0] * cosa) * vec2 + (-a_g[0,0] * cosa * sind - a_g[1,0] * sina * sind + a_g[2,0] * cosd) * vec3
   v = (a_g[0,1] * cosa * cosd + a_g[1,1] * sina * cosd + a_g[2,1] * sind) * vec1 + (-a_g[0,1] * sina + a_g[1,1] * cosa) * vec2 + (-a_g[0,1] * cosa * sind - a_g[1,1] * sina * sind + a_g[2,1] * cosd) * vec3
   w = (a_g[0,2] * cosa * cosd + a_g[1,2] * sina * cosd + a_g[2,2] * sind) * vec1 + (-a_g[0,2] * sina + a_g[1,2] * cosa) * vec2 + (-a_g[0,2] * cosa * sind - a_g[1,2] * sina * sind + a_g[2,2] * cosd) * vec3

   lsr_vel = numpy.array([-8.5, 13.38, 6.49]) # notice the sign of the first velocity 
   											#component, it is negative because 
									# in this program U points toward anticenter
   if (lsr is not None):   
      u = u + lsr_vel[0]
      v = v + lsr_vel[1]
      w = w + lsr_vel[2]
   
   return (u,v,w)
