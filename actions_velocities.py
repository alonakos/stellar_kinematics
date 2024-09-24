import sys
import glob
import math
import agama
import numpy as np
from astropy import units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from pygaia.astrometry.coordinates import CoordinateTransformation, Transformations
#from pygaia.astrometry.vectorastrometry import cartesianToSpherical, astrometryToPhaseSpace,phaseSpaceToAstrometry
from pygaia.astrometry.vectorastrometry import astrometryToPhaseSpace

#Set the units
agama.setUnits(mass=1,length=1,velocity=1)

#Parameters for debugging
#input_file_name = '../data/AEGIS/Original/gaia/aegis_rv/Trimmed_original_data_inputs_rv_aegis_gaia_dist_1.csv'
#potential = agama.Potential("../potentials/McMillan17.ini")

#Set the input file name
# input_file_name = glob.glob('temporary*')[0]

#Set the output file name
output_file_name = input_file_name[:-4] + '_outputs.csv'

#Get the potential and set the action finder
potential = agama.Potential("../../potentials/McMillan17.ini")
act_finder = agama.ActionFinder(potential)

def ICRS_to_GAL(phi, theta, muphistar, mutheta):
        """
        phi       - The longitude-like angle of the position of the source (radians)
        theta     - The latitude-like angle of the position of the source (radians)
        muphistar - Value of the proper motion in the longitude-like angle, multiplied by cos(latitude).
        mutheta   - Value of the proper motion in the latitude-like angle.
        """
    # transformation ICRS to GAL
    ctICRS2GAL =CoordinateTransformation(Transformations.ICRS2GAL)
        #
        ell,         bee = ctICRS2GAL.transformSkyCoordinates(phi,theta)
        muellstar, mubee = ctICRS2GAL.transformProperMotions (phi,theta, muphistar,mutheta)

        return ell,bee,muellstar,mubee

def astrometric_2_cartesian(DM, ell_radian, bee_radian, HRV, muellstar, mubee):
    parallax_mas = np.power(10., 2.-DM/5.)
    xHelio_pc, yHelio_pc, zHelio_pc, vxHelio, vyHelio, vzHelio = astrometryToPhaseSpace(ell_radian, bee_radian, parallax_mas, muellstar, mubee, HRV)
    xHelio_kpc = xHelio_pc/1000.
    yHelio_kpc = yHelio_pc/1000.
    zHelio_kpc = zHelio_pc/1000.

    return xHelio_kpc, yHelio_kpc, zHelio_kpc, vxHelio, vyHelio, vzHelio

def cartesian_2_cylindrical(x,y,z,vx,vy,vz):
    R =np.sqrt(x**2 + y**2)
    r =np.sqrt(x**2 + y**2 + z**2)

    phi=np.arctan2(y,x)
    vR     = (x*vx + y*vy)/R
    vTHETA = (R*vR - z*vz)/r
    vPHI   = (x*vy - y*vx)/R

    return R, phi, vR, vTHETA, vPHI

#calculate orbital parameters
def calculate_orbital_parameters(xGC_norm, yGC_norm, zGC_norm, vxGC_norm, vyGC_norm, vzGC_norm, inttime, numsteps):
        # galpy convention
        x_agama  =   xGC_norm
        y_agama  =   yGC_norm
        z_agama  =   zGC_norm
        vx_agama =   vxGC_norm
        vy_agama =   vyGC_norm
        vz_agama =   vzGC_norm

        #R_agama   =  np.sqrt(x_agama**2 + y_agama**2)
        #vR_agama  =  ( x_agama*vx_agama + y_agama*vy_agama )/R_agama
        #vT_agama  = -( x_agama*vy_agama - y_agama*vx_agama )/R_agama
        #phi_agama =  np.arctan2(y_agama,x_agama)


        #inttime=100.
        #numsteps=3000
    times = np.linspace(0, inttime, numsteps)
        times_c, c_orb_car = agama.orbit(ic=[x_agama,y_agama,z_agama,vx_agama,vy_agama,vz_agama], potential=potential, time=inttime, trajsize=numsteps)

    x = c_orb_car[:,0]
        y = c_orb_car[:,1]
        z = c_orb_car[:,2]
        vx= c_orb_car[:,3]
        vy= c_orb_car[:,4]
        vz= c_orb_car[:,5]
        R = np.sqrt(x**2 + y**2)
        r = np.sqrt(x**2 + y**2 + z**2)

        rmin = np.min(r)
        rmax = np.max(r)
        zmax = np.max(np.fabs(z))
        ecc  = (rmax-rmin)/(rmax+rmin)
        return rmin, rmax, zmax, ecc

#position and velocity of the Sun
_xGC_sun_kpc = -8.249
_yGC_sun_kpc =  0.
_zGC_sun_kpc =  0.
_vxGC_sun    =  11.1
_vyGC_sun    =  250.70
_vzGC_sun    =  7.25


#Open the data file
data = []
with open(input_file_name,'r') as file:
    data = file.readlines()

#Obtain the RA, DEC, distance, radial velocity, pmRA, and pmDEC from the file
#name list contains the stars which can be run through AGAMA while name_starred list are the stars
#which have a '*' in the data indicating that they cannot have orbits calculated due to missing information
#The name_all list is a list of all the stars that are found in the input file
name = []
name_all = []
name_starred = []
ra_dec = []
dist = []
dist_error = []
rv = []
rv_error = []
pmra = []
pmra_error = []
pmdec = []
pmdec_error = []
pmra_pmdec_corr_error = []
counter = -1
for line in data:
    counter = counter + 1
    if counter != 0:
        split_line = line.split(',')
        name_all.append(split_line[0])
        if "*" in split_line:
            name_starred.append(split_line[0])
        else:
            name.append(split_line[0])
            ra_dec.append(split_line[1]+' '+split_line[2])
            if float(split_line[3]) < 0:
                print split_line[0]+' has a distance less than zero!'+split_line[3]
            dist.append(float(split_line[3]))
            dist_error.append(float(split_line[4]))
            rv.append(float(split_line[5]))
            rv_error.append(float(split_line[6]))
            pmra.append(float(split_line[7]))
            pmra_error.append(float(split_line[8]))
            pmdec.append(float(split_line[9]))
            pmdec_error.append(float(split_line[10]))
            pmra_pmdec_corr_error.append(float(split_line[11]))

#print name, ra_dec, dist, dist_error, rv, rv_error, pmra, pmra_error, pmdec, pmdec_error, pmra_pmdec_corr_error

#Put the coordinates into a bin to get them easier in degrees and radians if there are any
ra_rad = []
dec_rad = []
if name:
    coordinates = SkyCoord(ra_dec,unit=(u.hourangle,u.deg))
    for coord in coordinates:
        ra_rad.append(coord.ra.rad)
        dec_rad.append(coord.dec.rad)

###############################################################################
#Calculate the actions and velocities and their errors using monte carlo methods

#Set the length and of the simulation
num_stars = len(name)
_N_MonteCarlo = 1000
np.random.seed()  #123 for debug purposes

#Arrays to put the results of the montecarlo simulation into
JrJzJphiE_Nstar_Nmontecarlo = np.ones((4, num_stars, _N_MonteCarlo))
XYZVxVy_Nstar_Nmontecarlo = np.ones((5, num_stars, _N_MonteCarlo))
vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo = np.ones((8, num_stars, _N_MonteCarlo))

Jr_mean = []
Jr_std = []
Jz_mean = []
Jz_std = []
Jphi_mean = []
Jphi_std = []
Jx_mean = []
Jx_std = []
Jy_mean = []
Jy_std = []
vr_mean = []
vr_std = []
vx_mean = []
vx_std = []
vy_mean = []
vy_std = []
vz_mean = []
vz_std = []
vphi_mean = []
vphi_std = []
Lx_mean = []
Lx_std = []
Ly_mean = []
Ly_std = []
energy_mean = []
energy_std = []
phi_mean = []
phi_std = []
rperi_mean = []
rperi_std = []
rapo_mean = []
rapo_std = []
zmax_mean = []
zmax_std = []
ecc_mean = []
ecc_std = []
for i in np.arange(0,num_stars):
    for j in np.arange(0,_N_MonteCarlo):
        RA_rad_ij = ra_rad[i] #no error
        DEC_rad_ij = dec_rad[i] #no error
        HRV_kms_ij = rv[i] + np.random.randn()*rv_error[i] #put in random error
        mean_pmRA_pmDEC = [pmra[i],pmdec[i]]
        #Get the error for proper motion RA and DEC
        e_pmRA = pmra_error[i]
        e_pmDEC = pmdec_error[i]
        e_pmRA_pmDEC = pmra_pmdec_corr_error[i]
        #Calculate the sigma from proper motion RA and DEC
        Sigma_pmRA = e_pmRA**2
        Sigma_pmDEC = e_pmDEC**2
        Sigma_pmRA_pmDEC = e_pmRA_pmDEC*e_pmRA*e_pmDEC
        Sigma = [[Sigma_pmRA, Sigma_pmRA_pmDEC],
             [Sigma_pmRA_pmDEC, Sigma_pmDEC]]
        #Create an array of the proper motion RA and DEC
        pmRAstar_masyr_ij_array, pmDECstar_masyr_ij_array = np.random.multivariate_normal(mean=mean_pmRA_pmDEC, cov=Sigma, size=1).T
        pmRAstar_masyr_ij = pmRAstar_masyr_ij_array[0]
        pmDECstar_masyr_ij = pmDECstar_masyr_ij_array[0]
        #Create a distance with error
        dist_kpc_ij = dist[i] + np.random.rand()*dist_error[i]
        while dist_kpc_ij < 0.0:
            dist_kpc_ij = dist[i] + np.random.rand()*dist_error[i]
        parallax_mas_ij = 1.0/dist_kpc_ij
        #Calculate a magnitude from the calculated parallax
            DM_mag_ij = 5.*(np.log10(1./parallax_mas_ij) + 2.0)
        #Put the star into galactic coordinates
        l_rad_ij, b_rad_ij, pmlstar_masyr_ij, pmbstar_masyr_ij = ICRS_to_GAL(RA_rad_ij, DEC_rad_ij, pmRAstar_masyr_ij, pmDECstar_masyr_ij)
        #Convert to Cartesian
        xH, yH, zH, vxH, vyH, vzH = astrometric_2_cartesian(DM_mag_ij, l_rad_ij, b_rad_ij, HRV_kms_ij, pmlstar_masyr_ij, pmbstar_masyr_ij)
        #Normalize the position
        xGC_norm = (xH + _xGC_sun_kpc)
        yGC_norm = (yH + _yGC_sun_kpc)
        zGC_norm = (zH + _zGC_sun_kpc)
         #normalize  velocity
                vxGC_norm = (vxH + _vxGC_sun)
                vyGC_norm = (vyH + _vyGC_sun)
                vzGC_norm = (vzH + _vzGC_sun)
        #cylindrical coordinate
                R, phi, vR, vTHETA, vPHI = cartesian_2_cylindrical(xGC_norm, yGC_norm, zGC_norm, vxGC_norm, vyGC_norm, vzGC_norm)
                #vPERP = np.sqrt(vR**2 + vzGC_norm**2)
        #Calculate the energy
        energy = potential.potential([xGC_norm,yGC_norm,zGC_norm]) + 0.5*(vxGC_norm**2+vyGC_norm**2+vzGC_norm**2)
        #Calculate the actions
        xv6d = np.array( [xGC_norm, yGC_norm, zGC_norm, vxGC_norm, vyGC_norm, vzGC_norm]).T
        Jr, Jz, Jphi = act_finder(xv6d)

        #Since we are using a right-handed coordinate system
        #The phi results need to have sign switch in order
        #for prograde stars to have a positive action
        vPHI = -vPHI
        Jphi = -Jphi

        rmin, rmax, zmax, ecc= np.nan, np.nan, np.nan, np.nan
           if (np.isnan(Jr)):
                    Jr = -1.
                if (np.isnan(Jz)):
                    Jz = -1.
                if (np.isnan(energy)):
            energy = 1
        #Get the orbital parameters
                if (energy<0.):
                    rmin, rmax, zmax, ecc = calculate_orbital_parameters(xGC_norm,yGC_norm,zGC_norm,vxGC_norm,vyGC_norm,vzGC_norm,100,10001)
            print i, j
        else:
            Jr, Jz, Jphi, energy, vR, vzGC_norm, vPHI, phi = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            print "here"
        #Put the actions and the energy into an array
        JrJzJphiE_Nstar_Nmontecarlo[0,i,j] = Jr
		JrJzJphiE_Nstar_Nmontecarlo[1,i,j] = Jz
		JrJzJphiE_Nstar_Nmontecarlo[2,i,j] = Jphi
		JrJzJphiE_Nstar_Nmontecarlo[3,i,j] = energy
        #Put the velocities into an array
        vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i,j] = vR
		vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i,j] = vzGC_norm
		vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[2,i,j] = vPHI
		vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i,j] = phi
        #Put the positions and remaining velocities into an array
        XYZVxVy_Nstar_Nmontecarlo[0,i,j] = xGC_norm
        XYZVxVy_Nstar_Nmontecarlo[1,i,j] = yGC_norm
        XYZVxVy_Nstar_Nmontecarlo[2,i,j] = zGC_norm
        XYZVxVy_Nstar_Nmontecarlo[3,i,j] = vxGC_norm
        XYZVxVy_Nstar_Nmontecarlo[4,i,j] = vyGC_norm
        #Put the distances and eccentricities into an array
		vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[4,i,j] = rmin
		vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[5,i,j] = rmax
		vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[6,i,j] = zmax
		vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[7,i,j] = ecc

        #Flush the output to the output file
        sys.stdout.flush()

    #Calculate the mean and std of the actions
    Jr_mean.append(np.nanmean(JrJzJphiE_Nstar_Nmontecarlo[0,i]))
    Jr_std.append(np.nanstd(JrJzJphiE_Nstar_Nmontecarlo[0,i]))
    Jz_mean.append(np.nanmean(JrJzJphiE_Nstar_Nmontecarlo[1,i]))
        Jz_std.append(np.nanstd(JrJzJphiE_Nstar_Nmontecarlo[1,i]))
    Jphi_mean.append(np.nanmean(JrJzJphiE_Nstar_Nmontecarlo[2,i]))
        Jphi_std.append(np.nanstd(JrJzJphiE_Nstar_Nmontecarlo[2,i]))
        #Calculate the mean and std of the energy
    energy_mean.append(np.nanmean(JrJzJphiE_Nstar_Nmontecarlo[3,i]))
        energy_std.append(np.nanstd(JrJzJphiE_Nstar_Nmontecarlo[3,i]))
    #Calculate the mean and std of the velocities
    vr_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i]))
    vr_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i]))
    vz_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i]))
    vz_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i]))
    vphi_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[2,i]))
    vphi_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[2,i]))
    #Calculate the mean and std of the angular momenta
    Lx_mean.append(np.nanmean(XYZVxVy_Nstar_Nmontecarlo[1,i]*vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i] - XYZVxVy_Nstar_Nmontecarlo[2,i]*XYZVxVy_Nstar_Nmontecarlo[4,i]))
    Lx_std.append(np.nanstd(XYZVxVy_Nstar_Nmontecarlo[1,i]*vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i] - XYZVxVy_Nstar_Nmontecarlo[2,i]*XYZVxVy_Nstar_Nmontecarlo[4,i]))
    Ly_mean.append(np.nanmean(XYZVxVy_Nstar_Nmontecarlo[2,i]*XYZVxVy_Nstar_Nmontecarlo[3,i] - XYZVxVy_Nstar_Nmontecarlo[0,i]*vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i]))
    Ly_std.append(np.nanstd(XYZVxVy_Nstar_Nmontecarlo[2,i]*XYZVxVy_Nstar_Nmontecarlo[3,i] - XYZVxVy_Nstar_Nmontecarlo[0,i]*vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i]))
    #Calculate the mean and std of the distances and eccentricities
    rperi_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[4,i]))
    rperi_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[4,i]))
    rapo_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[5,i]))
    rapo_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[5,i]))
    zmax_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[6,i]))
    zmax_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[6,i]))
    ecc_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[7,i]))
    ecc_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[7,i]))
    #Calculate the mean and std of phi (Used to convert to Cartesian Coordinates)
    phi_mean.append(np.nanmean(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i]))
    phi_std.append(np.nanstd(vr_vz_vphi_phi_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i]))
    #Calculate the mean and std of Cartesian actions and velocities
    Jx_mean.append(np.nanmean(Jr_mean[i]*np.cos(phi_mean[i])))
    Jx_std.append(np.sqrt(np.cos(phi_mean[i])**2*Jr_std[i]*2+Jr_mean[i]**2*np.sin(phi_mean[i])**2*phi_std[i]**2))
    Jy_mean.append(np.nanmean(Jr_mean[i]*np.sin(phi_mean[i])))
    Jy_std.append(np.sqrt(np.sin(phi_mean[i])**2*Jr_std[i]**2+Jr_mean[i]**2*np.cos(phi_mean[i])**2*phi_std[i]**2))
    vx_mean.append(np.nanmean(XYZVxVy_Nstar_Nmontecarlo[3,i]))
    vx_std.append(np.nanstd(XYZVxVy_Nstar_Nmontecarlo[3,i]))
    vy_mean.append(np.nanmean(XYZVxVy_Nstar_Nmontecarlo[4,i]))
    vy_std.append(np.nanstd(XYZVxVy_Nstar_Nmontecarlo[4,i]))

#Put the velocities, actions, energy, Apocentric/Pericentric Distance, z max and the eccentricity into output file
with open(output_file_name,'w') as file:
    file.write('Name,v_r(km/s),v_r_std(km/s),v_phi(km/s),v_phi_std(km/s),v_x(km/s),v_x_std(km/s),v_y(km/s),v_y_std(km/s),v_z(km/s),v_z_std(km/s),j_r(kpc km s^-1),j_r_std(kpc km s^-1),j_phi(kpc s^(-1)),j_phi_std(kpc km s^(-1)),j_x(kpc km s^(-1)),j_x_std(kpc km s^(-1)),j_y(kpc km s^(-1)),j_y_std(kpc km s^(-1)),j_z(kpc km s^(-1)),j_z_std(kpc km s^(-1)),L_x(km s^(-1) kpc),L_x_std(km s^(-1) kpc),L_y(km s^(-1) kpc),L_y_std(km s^(-1) kpc),energy_mean(km^2 s^(-2)),energy_std(km^2 s^(-2)),rperi_mean(kpc),rperi_std(kpc),rapo_mean(kpc),rapo_std(kpc),zmax_mean(kpc),zmax_std(kpc),ecc_mean,ecc_std\n')
    for i in np.arange(0,len(name_all)):
        #Check to see if AGAMA ran for the star
        if name_all[i] in name and not name_all[i] in name_starred:
            #Get the index of the star in the good name list
            name_index = name.index(name_all[i])
            #Write the name of the object
            file.write(name[name_index]+',')
            #Write the velocities
            file.write(str(vr_mean[name_index])+',')
            file.write(str(vr_std[name_index])+',')
            file.write(str(vphi_mean[name_index])+',')
            file.write(str(vphi_std[name_index])+',')
            file.write(str(vx_mean[name_index])+',')
            file.write(str(vx_std[name_index])+',')
            file.write(str(vy_mean[name_index])+',')
            file.write(str(vy_std[name_index])+',')
            file.write(str(vz_mean[name_index])+',')
            file.write(str(vz_std[name_index])+',')
            #Write the actions
            file.write(str(Jr_mean[name_index])+',')
            file.write(str(Jr_std[name_index])+',')
            file.write(str(Jphi_mean[name_index])+',')
            file.write(str(Jphi_std[name_index])+',')
            file.write(str(Jx_mean[name_index])+',')
            file.write(str(Jx_std[name_index])+',')
            file.write(str(Jy_mean[name_index])+',')
            file.write(str(Jy_std[name_index])+',')
            file.write(str(Jz_mean[name_index])+',')
            file.write(str(Jz_std[name_index])+',')
            #Write the angular momenta
            file.write(str(Lx_mean[name_index])+',')
            file.write(str(Lx_std[name_index])+',')
            file.write(str(Ly_mean[name_index])+',')
            file.write(str(Ly_std[name_index])+',')
            #Write the energy
            file.write(str(energy_mean[name_index])+',')
            file.write(str(energy_std[name_index])+',')
            #Write the Apocentric/Pericentric Distance, z max and the eccentricity
            file.write(str(rperi_mean[name_index])+',')
            file.write(str(rperi_std[name_index])+',')
            file.write(str(rapo_mean[name_index])+',')
            file.write(str(rapo_std[name_index])+',')
            file.write(str(zmax_mean[name_index])+',')
            file.write(str(zmax_std[name_index])+',')
            file.write(str(ecc_mean[name_index])+',')
            file.write(str(ecc_std[name_index])+'\n')
        #Check to see if there were no parameters for AGAMA to run
        elif name_all[i] in name_starred and not name_all[i] in name:
            file.write(name_all[i]+',*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*\n')
        #Make sure nothing went wrong
        else:
            file.write(name_all[i]+',SOMETHING WENT WRONG! STAR NOT FOUND!')
