import sys
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
#sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
sys.path.append('/home/kkosciw/Python_Scripts/arepo_package')
import site
site.addsitedir("/home/kkosciw/Python_Scripts")


import arepo_package
#%pylab inline

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter,gaussian_filter1d

import astro_constants as ac



def find_BH_position(basePath,desired_redshift,p_id,choose_most_massive=0): 
    bh_IDs,output_redshift=arepo_package.get_particle_property(basePath,'ParticleIDs',5,desired_redshift,list_all=False)
    bh_positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',5,desired_redshift,list_all=False)
 
    bh_masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',5,desired_redshift,list_all=False)
    
    if (choose_most_massive):
        center=bh_positions[bh_masses==numpy.amax(bh_masses)]
    else:
        center=(bh_positions[bh_IDs==p_id])
    return center

def find_FOF_position(basePath,desired_redshift,choose_most_massive=0):
    fofpos,output_redshift=arepo_package.get_group_property(basePath,'GroupPos',desired_redshift,list_all=False)    
    fofmass,output_redshift=arepo_package.get_group_property(basePath,'GroupMass',desired_redshift,list_all=False)
    if (choose_most_massive):
        maxpos=numpy.where(fofmass==numpy.amax(fofmass))
        center=fofpos[maxpos,:]
    return center

def find_SubHalo_position(basePath,desired_redshift,desired_index):
    SubhaloGrNr,output_redshift=arepo_package.get_subhalo_property(basePath,'SubhaloGrNr',desired_redshift,list_all=False)
    SubhaloPos,output_redshift=arepo_package.get_subhalo_property(basePath,'SubhaloPos',desired_redshift,list_all=False)
    Subhalo_Indices=numpy.arange(len(SubhaloGrNr))
    SH_Index=Subhalo_Indices[desired_index]
    center=SubhaloPos[SH_Index]
    center=numpy.reshape(center,(1,3))
    #Can use SubhaloVmaxRad for analysis purposes
    SubhaloVmaxRad,output_redshift=arepo_package.get_subhalo_property(basePath,'SubhaloVmaxRad',desired_redshift,list_all=False)
    SubhaloRad=SubhaloVmaxRad[SH_Index]

    
    return center,SubhaloRad

def orient_plane(positions,perpendicular_vector):
    #Orienting the galaxy to be face-on or edge-on
    def mag(vector):
        return numpy.sqrt(sum(vector**2))    
    A=numpy.array([1,1,1])
    unit_A=A/mag(A)
    unit_perpendicular_vector=perpendicular_vector/mag(perpendicular_vector)
    
    L=unit_perpendicular_vector      
    LxA=numpy.cross(L,unit_A)
    LxLxA=numpy.cross(L,LxA)
    
    LxA/=mag(LxA)
    LxLxA/=mag(LxLxA)
    
    zpos_new=numpy.dot(positions,L)
    
    xpos_new=numpy.dot(positions,LxA)
    
    ypos_new=numpy.dot(positions,LxLxA)

    new_positions=numpy.transpose(numpy.array([xpos_new,ypos_new,zpos_new]))
    
    return new_positions

def extract_slice(basePath,p_type,desired_center,desired_redshift,field,planeofsky_dimensions=(100,100),lineofsight_dimension=20,plane='xy',orient=0):
    #Gets information on which particle property is desired and masks out any particles that are not within the box
    positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',p_type,desired_redshift,list_all=False)       
    #Subhalo Positions
    #positions,output_redshift=arepo_package.get_subhalo_property(basePath,'SubhaloPos',desired_redshift,list_all=False)    
    print(positions.shape)
    if (field=='SHSFR'):
        SubhaloSFR,output_redshift=arepo_package.get_subhalo_property(basePath,'SubhaloSFR',desired_redshift,list_all=False)
        particle_property=SubhaloSFR
    if (field=='Density'):
        masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
        masses*=1e10 
        particle_property=masses
    if (field=='Metallicity'):
        metallicity,output_redshift=arepo_package.get_particle_property(basePath,'GFM_Metallicity',p_type,desired_redshift,list_all=False)  
        metallicity/=0.0127
        particle_property=metallicity        
    if (field=='Temperature'):
        internal_energy,output_redshift=arepo_package.get_particle_property(basePath,'InternalEnergy',p_type,desired_redshift,list_all=False)
        electron_abundance,output_redshift=arepo_package.get_particle_property(basePath,'ElectronAbundance',p_type,desired_redshift,list_all=False)
        g_minus_1 = (5.0/3.0) - 1.0
        XH = 0.76
        mu=(4*ac.MP)/(1+3*XH+4*XH*electron_abundance)
        gas_temperature = g_minus_1*(internal_energy/ac.KB)*(10**10)*mu
        particle_property=gas_temperature
    if (field=='Velocity Magnitude'):
        velocity,output_redshift=arepo_package.get_particle_property(basePath,'Velocities',p_type,desired_redshift,list_all=False) 
        xvel=velocity[:,0]
        yvel=velocity[:,1]
        zvel=velocity[:,2]
        mag_vel=numpy.sqrt((xvel**2)+(yvel**2)+(zvel**2))
        particle_property=mag_vel
    print(output_redshift)
    positions_relative_to_center=positions-numpy.ravel(desired_center)
      
    if (orient):
        positions_relative_to_center=orient_plane(positions_relative_to_center,perpendicular_vector)
         
    if (plane=='xy'):
        planeofsky_pos1=positions_relative_to_center[:,0]
        planeofsky_pos2=positions_relative_to_center[:,1]
        lineofsight_pos=positions_relative_to_center[:,2]
        
    if (plane=='xz'):
        planeofsky_pos1=positions_relative_to_center[:,0]
        planeofsky_pos2=positions_relative_to_center[:,2]
        lineofsight_pos=positions_relative_to_center[:,1]
    
    if (plane=='yz'):
        planeofsky_pos1=positions_relative_to_center[:,1]
        planeofsky_pos2=positions_relative_to_center[:,2]
        lineofsight_pos=positions_relative_to_center[:,0]
        
    observable_positions=numpy.transpose(numpy.array([planeofsky_pos1,planeofsky_pos2,lineofsight_pos]))
    
    mask_planeofsky_1=numpy.abs(planeofsky_pos1)<planeofsky_dimensions[0]/2
    mask_planeofsky_2=numpy.abs(planeofsky_pos2)<planeofsky_dimensions[1]/2
    mask_lineofsight=numpy.abs(lineofsight_pos)<lineofsight_dimension/2
    mask=(mask_planeofsky_1&mask_planeofsky_2)&mask_lineofsight
        
    return observable_positions[mask],particle_property[mask],output_redshift

def visualize(final_positions,final_property,number_of_pixels,field):
    n_planeofsky1=number_of_pixels
    n_planeofsky2=number_of_pixels
    planeofsky_pos1=final_positions[:,0]
    planeofsky_pos2=final_positions[:,1]
    lineofsight_pos=final_positions[:,2]
    

    
    planeofsky1_grid=numpy.linspace(min(planeofsky_pos1),max(planeofsky_pos1),n_planeofsky1)
    planeofsky2_grid=numpy.linspace(min(planeofsky_pos2),max(planeofsky_pos2),n_planeofsky2)
    pixelsize_planeofsky1=numpy.diff(planeofsky1_grid)[0]
    pixelsize_planeofsky2=numpy.diff(planeofsky2_grid)[0]
    pixelsize_lineofsight=numpy.amax(lineofsight_pos)-numpy.amin(lineofsight_pos)

    proj_property=[]
    if (field=='Density'):
        pixel_volume=pixelsize_lineofsight*pixelsize_planeofsky1*pixelsize_planeofsky2
 
    First,Second=numpy.meshgrid(planeofsky1_grid,planeofsky2_grid)
    for pixel_position2 in planeofsky2_grid:
        for pixel_position1 in planeofsky1_grid:
            mask1=(planeofsky_pos1>(pixel_position1-pixelsize_planeofsky1/2.0))&(planeofsky_pos1<(pixel_position1+pixelsize_planeofsky1/2.))
            mask2=(planeofsky_pos2>(pixel_position2-pixelsize_planeofsky2/2.0))&(planeofsky_pos2<(pixel_position2+pixelsize_planeofsky2/2.))
       
            mask=(mask1) & (mask2)
            if (field=='Density'):
                proj_property.append(numpy.sum(final_property[mask])/pixel_volume)
            else:
                proj_property.append(numpy.average(final_property[mask]))
    #The Density is the only property that won't have nans, the 1e-19 will not be plotted later, just useful for when the log needs to be taken                      
    proj_property=numpy.asarray(proj_property)
    proj_property[numpy.isnan(proj_property)]=1e-19
    Proj_property=proj_property.reshape(number_of_pixels,number_of_pixels)
    #Helps with making the plots look better. Higher sigma = more smoothed plot
    Proj_property=gaussian_filter(Proj_property,sigma=1)
       

    if (field=='Density'):
        print("making density")
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(proj_property[proj_property>0]),vmax=Proj_property.max()),cmap='Greys_r')
        #Allows for multiple plots to be on same scale
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=1e3,vmax=1e7),cmap='Greys_r')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas density ($M_{\odot}/(cKpc/h)^3$)',fontsize=15)
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)
        #bh=plt.Circle((bhcoord[0],bhcoord[1]),0.5,color='black')
        #ax.add_artist(bh)
  
    if (field=='Metallicity'):
        print('making metal')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='plasma')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=1e-8,vmax=1e0),cmap='plasma')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas metallicity $(Fe/H)/(Fe/H)_{\odot}$',fontsize=15)
        #cbar.ax.set_ylim([cbar.norm(10e1),cbar.norm(10e-7)])
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)

    if (field=='Temperature'):
        print('making temperature')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='hot')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=1e2,vmax=1e7),cmap='hot')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas temperature $K$',fontsize=15)
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)
 
    if (field=='Velocity Magnitude'):
        print('making velocity mag')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='viridis')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=1e1,vmax=1e3),cmap='viridis')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas velocity magnitude $km \sqrt{a}/s$',fontsize=15)
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)

    
def stellar_ang_mom(basePath,desired_redshift,desired_center,box_length=40):
    #Calculating the angular momentum allows for the galaxy to be oriented face-on or edge-on
    p_type=4
    positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',p_type,desired_redshift,list_all=False)       
    masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
    velocities,output_redshift=arepo_package.get_particle_property(basePath,'Velocities',p_type,desired_redshift,list_all=False)
    masses*=1e10
    velocities/=3.0886e16
    
    xvel=velocities[:,0] #units switched from km/s to kpc/s
    yvel=velocities[:,1]
    zvel=velocities[:,2]

    xpos=positions[:,0]
    ypos=positions[:,1]
    zpos=positions[:,2]
    
    center=desired_center[0]
    center = numpy.ravel(center) 
    mask_x=numpy.abs(xpos-center[0])<box_length/2
    mask_y=numpy.abs(ypos-center[1])<box_length/2
    mask_z=numpy.abs(zpos-center[2])<box_length/2
    mask=(mask_x&mask_y)&mask_z
    xvel=xvel[mask]
    yvel=yvel[mask]
    zvel=zvel[mask]
    
    xpos=xpos[mask]
    ypos=ypos[mask]
    zpos=zpos[mask]

    masses=masses[mask]
    
    #Calculate angular momentum
    Lx=(ypos*zvel-zpos*yvel)*masses
    Lx=sum(Lx)
    Ly=(xpos*zvel-zpos*xvel)*masses
    Ly=sum(Ly)
    Lz=(xpos*yvel-ypos*xvel)*masses
    Lz=sum(Lz)    
    L_total=numpy.transpose(numpy.array([Lx,Ly,Lz]))
     
    p_particles=numpy.transpose(numpy.array([xvel*masses,yvel*masses,zvel*masses]))


    p_gal=numpy.sum(p_particles,axis=0)

    num=numpy.transpose(numpy.array([xpos*masses,ypos*masses,zpos*masses]))
    com=numpy.sum(num,axis=0)/numpy.sum(masses)
    

    L_orbital=numpy.cross(com,p_gal)

    L_spin=L_total-L_orbital


    return L_spin

def line_plots(basePath,desired_redshift,p_type,field):
    if (field=='Metallicity'):
        bh_IDs,output_redshift=arepo_package.get_particle_property(basePath,'ParticleIDs',5,desired_redshift,list_all=False)
        desired_center=find_BH_position(basePath,desired_redshift,p_id=bh_IDs[BH_index],choose_most_massive=0)
        pos,metallicity,output_redshift=extract_slice(basePath,p_type,desired_center,desired_redshift,field,planeofsky_dimensions=(100,100),lineofsight_dimension=20,plane='xy',orient=0)

    
    return metallicity,output_redshift
 
#Old Sims
#path_to_uniform_run='/blue/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
#New Sims
path_to_uniform_run='/blue/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_ZOOM_RUNS4/'

#Old Sims
#uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
#uniform_run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
#uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass5.90_logFOFseedmass10.70/AREPO'
#basePath=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
#NewSims
#uniform_run='density_and_metallicity_based_criterion_zoom_levelmin7_levelmax9_haloindex4_redshift5.00_logbhseedmass5.90_NSC/AREPO/'
#uniform_run='density_and_metallicity_based_criterion_zoom_levelmin7_levelmax10_haloindex4_redshift5.00_logbhseedmass5.00_NSC/AREPO/'
uniform_run='density_and_metallicity_based_criterion_zoom_levelmin7_levelmax11_haloindex4_redshift5.00_logbhseedmass4.10_NSC/AREPO/'
#basePath=path_to_uniform_run+uniform_run+'output_upto_4'
#Use for Lev Max 11 only
basePath=path_to_uniform_run+uniform_run+'output_upto_4_previous_version'



#Make the profile plots

#desired_redshift=5
#BH_index=0
#desired_index=0
plane='xy'

#bh_IDs,output_redshift=arepo_package.get_particle_property(basePath,'ParticleIDs',5,desired_redshift,list_all=False)
#desired_center=find_BH_position(basePath,desired_redshift,p_id=bh_IDs[BH_index],choose_most_massive=1)
#desired_center=find_FOF_position(basePath,desired_redshift,choose_most_massive=1)
#desired_center,shrad=find_SubHalo_position(basePath,desired_redshift,desired_index)



desired_redshift=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
desired_index=[0,1,0,328,0,0,0,0,123,38,0,13,16,1,1,1,2,7,10]
desired_center,maxrad=find_SubHalo_position(basePath,5,0)
halfrad=maxrad/2
qrtrad=maxrad/4

final_prop_full=[]
final_prop_half=[]
final_prop_qrt=[]
for  dr,shind in zip(desired_redshift,desired_index):
    desired_center,shrad=find_SubHalo_position(basePath,dr,shind)
    perpendicular_vector=stellar_ang_mom(basePath,dr,desired_center)
    final_positions,final_property,output_redshift=extract_slice(basePath,0,desired_center,dr,'Density',planeofsky_dimensions=(200,200),lineofsight_dimension=200,plane=plane,orient=1)

    #print(final_positions.shape)
    part_rad=numpy.sqrt((final_positions[:,0]**2)+(final_positions[:,1]**2)+(final_positions[:,2]**2))
    rad_mask_full=part_rad<maxrad

    xposrad = final_positions[:,0]
    yposrad = final_positions[:,1]
    zposrad = final_positions[:,2]

    xposrad = xposrad[rad_mask_full]
    yposrad = yposrad[rad_mask_full]
    zposrad = zposrad[rad_mask_full]

    #final_prop_full.append(numpy.average(final_property[rad_mask_full]))
    #Use for Density only
    final_prop_full.append(numpy.average(numpy.sum(final_property[rad_mask_full])/(4/3*numpy.pi*(maxrad**3))))
   
    rad_mask_half=part_rad<halfrad

    xposrad = final_positions[:,0]
    yposrad = final_positions[:,1]
    zposrad = final_positions[:,2]

    xposrad = xposrad[rad_mask_half]
    yposrad = yposrad[rad_mask_half]
    zposrad = zposrad[rad_mask_half]

    #final_prop_half.append(numpy.average(final_property[rad_mask_half]))
    #Use for Density only
    final_prop_half.append(numpy.average(numpy.sum(final_property[rad_mask_half])/(4/3*numpy.pi*(halfrad**3))))

    rad_mask_qrt=part_rad<qrtrad

    xposrad = final_positions[:,0]
    yposrad = final_positions[:,1]
    zposrad = final_positions[:,2]

    xposrad = xposrad[rad_mask_qrt]
    yposrad = yposrad[rad_mask_qrt]
    zposrad = zposrad[rad_mask_qrt]

    #final_prop_qrt.append(numpy.average(final_property[rad_mask_qrt]))
    #Use for Density only
    final_prop_qrt.append(numpy.average(numpy.sum(final_property[rad_mask_qrt])/(4/3*numpy.pi*(qrtrad**3))))




fig, ax = plt.subplots(1,1)
#ax.set_facecolor('xkcd:black')
#visualize(final_positions,final_property,300,'Density')
#ax.tick_params(labelsize=20)
ax.set_ylabel('log Density $M_{\odot}/(cKpc/h)^3$')
ax.set_yscale('log')
ax.plot(desired_redshift,final_prop_full,label='Radius=%.1f'%maxrad,color='red')
ax.plot(desired_redshift,final_prop_half,label='Radius=%.1f'%halfrad,color='blue')
ax.plot(desired_redshift,final_prop_qrt,label='Radius=%.1f'%qrtrad,color='green')
ax.set_xlabel('redshift')
plt.legend(loc='upper right')
plt.title('Avg Density Lev Max 11',size=15)
plt.tight_layout()
fig.savefig('Prof_Package_Test')
plt.close()



