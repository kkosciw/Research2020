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
    #print(bh_positions,bh_masses)
    if (choose_most_massive):
        center=bh_positions[bh_masses==numpy.amax(bh_masses)]
    else:
        center=(bh_positions[bh_IDs==p_id])
    return center

def find_FOF_position(basePath,desired_redshift,choose_most_massive=0):
    fofpos,output_redshift=arepo_package.get_group_property(basePath,'GroupPos',desired_redshift,list_all=True)
    print(output_redshift)
    fofmass,output_redshift=arepo_package.get_group_property(basePath,'GroupMass',desired_redshift,list_all=False)
    if (choose_most_massive):
        maxpos=numpy.where(fofmass==numpy.amax(fofmass))
        center=fofpos[maxpos,:]
    return center

def orient_plane(positions,perpendicular_vector):
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
    print(new_positions)
    print('is this a metallica?')
    return new_positions

def extract_slice(basePath,p_type,desired_center,desired_redshift,field,planeofsky_dimensions=(1000,1000),lineofsight_dimension=20,plane='xy',orient=0):
    positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',p_type,desired_redshift,list_all=False)    
    if (field=='Density'):
        masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
        masses*=1e10 
        particle_property=masses
    if (field=='Metallicity'):
        metallicity,output_redshift=arepo_package.get_particle_property(basePath,'GFM_Metallicity',p_type,desired_redshift,list_all=False)  
        metallicity/=0.0127
        print(metallicity)
        particle_property=metallicity        
    if (field=='Temperature'):
        internal_energy,output_redshift=arepo_package.get_particle_property(basePath,'InternalEnergy',p_type,desired_redshift,list_all=False)
        electron_abundance,output_redshift=arepo_package.get_particle_property(basePath,'ElectronAbundance',p_type,desired_redshift,list_all=False)
        g_minus_1 = (5.0/3.0) - 1.0
        XH = 0.76
        mu=(4*ac.MP)/(1+3*XH+4*XH*electron_abundance)
        gas_temperature = g_minus_1*(internal_energy/ac.KB)*(10**10)*mu
        #print(gas_temperature)
        #print(output_redshift)
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
    #print(positions_relative_to_center[:,0])   
    if (orient):
        positions_relative_to_center=orient_plane(positions_relative_to_center,perpendicular_vector)
        #print('help me beech!') 
    if (plane=='xy'):
        planeofsky_pos1=positions_relative_to_center[:,0]
        #print(planeofsky_pos1)
        #print 'a'
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
    #print planeofsky_dimensions[0]
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
                #print(proj_property)                   
    proj_property=numpy.asarray(proj_property)
    print(proj_property)
    proj_property[numpy.isnan(proj_property)]=1e-19
    print(proj_property)
    Proj_property=proj_property.reshape(number_of_pixels,number_of_pixels)
    Proj_property=gaussian_filter(Proj_property,sigma=1)
       

    if (field=='Density'):
        print("making")
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(proj_property[proj_property>0]),vmax=Proj_property.max()),cmap='Greys_r')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=4e4,vmax=Proj_property.max()),cmap='Greys_r')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas density ($M_{\odot}/(cKpc/h)^3$)',fontsize=15)
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)
        #bh=plt.Circle((bhcoord[0],bhcoord[1]),0.5,color='black')
        #ax.add_artist(bh)
  
    if (field=='Metallicity'):
        print('making metal')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='plasma')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=1e-8,vmax=1e-1),cmap='plasma')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas metallicity $(Fe/H)/(Fe/H)_{\odot}$',fontsize=15)
        #cbar.ax.set_ylim([cbar.norm(10e1),cbar.norm(10e-7)])
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)

    if (field=='Temperature'):
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='hot')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=1e2,vmax=1e7),cmap='hot')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas temperature $K$',fontsize=15)
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)
 
    if (field=='Velocity Magnitude'):
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='viridis')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas velocity magnitude $km \sqrt{a}/s$',fontsize=15)
        ax.set_xlabel('Plane of sky 1',fontsize=15)
        ax.set_ylabel('Plane of sky 2',fontsize=15)

    
def stellar_ang_mom(basePath,desired_redshift,desired_center,box_length=40):
    p_type=4
    positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',p_type,desired_redshift,list_all=False)       
    masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
    #print('masses')
    #print(masses)
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
    #numpy.set_printoptions(threshold=numpy.inf)
    #print(numpy.abs(xpos-center[0]))
    #print(numpy.abs(ypos-center[1]))
    #print(numpy.abs(zpos-center[2]))
    #print(mask)
    xvel=xvel[mask]
    yvel=yvel[mask]
    zvel=zvel[mask]
    
    xpos=xpos[mask]
    ypos=ypos[mask]
    zpos=zpos[mask]

    masses=masses[mask]
    #print('mass mask',masses)    
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
    #print(num)
    #print('sup?')
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
 

#path_to_uniform_run='/blue/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
path_to_uniform_run='/blue/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_ZOOM_RUNS4/'
#uniform_run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
#uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass5.90_logFOFseedmass10.70/AREPO'

#uniform_run='density_and_metallicity_based_criterion_zoom_levelmin7_levelmax9_haloindex4_redshift5.00_logbhseedmass5.90_NSC/AREPO/'
#uniform_run='density_and_metallicity_based_criterion_zoom_levelmin7_levelmax10_haloindex4_redshift5.00_logbhseedmass5.00_NSC/AREPO/'
uniform_run='density_and_metallicity_based_criterion_zoom_levelmin7_levelmax11_haloindex4_redshift5.00_logbhseedmass4.10_NSC/AREPO/'

#basePath=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
#basePath=path_to_uniform_run+uniform_run+'output_upto_4'
basePath=path_to_uniform_run+uniform_run+'output_upto_4_previous_version'

desired_redshift=12
BH_index=0
#Lev Max 10 z 21 p_id=1001480254 ind=0
#Lev Max 10 z 20 p_id=1001538158 ind=1
#Lev Max 10 z 18 p_id=1004872384 ind=0
#p_id=1002600983
#p_id=1005831581
#p_id=1056102916
plane='xy'

bh_IDs,output_redshift=arepo_package.get_particle_property(basePath,'ParticleIDs',5,desired_redshift,list_all=False)
#print(bh_IDs)
desired_center=find_BH_position(basePath,desired_redshift,p_id=bh_IDs[BH_index],choose_most_massive=0)
print('dc',desired_center)
#desired_center=find_FOF_position(basePath,desired_redshift,choose_most_massive=1)
#print(desired_center)

perpendicular_vector=stellar_ang_mom(basePath,desired_redshift,desired_center)
final_positions,final_property,output_redshift=extract_slice(basePath,0,desired_center,desired_redshift,'Temperature',planeofsky_dimensions=(100,100),lineofsight_dimension=20,plane=plane,orient=1)
#print(final_positions)
fig, ax = plt.subplots(1,1,figsize=(11,9))

ax.set_facecolor('xkcd:black')
visualize(final_positions,final_property,300,'Temperature')
ax.tick_params(labelsize=20)
plt.title('Temeprature Level Max 11  z=%.1f'% output_redshift,size=15)
fig.savefig('Prof_Package_Test')


'''
red=0.2
fofmass,output_redshift=arepo_package.get_group_property(basePath,'GroupMass',red,list_all=True)
fofpos,output_redshift=arepo_package.get_group_property(basePath,'GroupPos',red,list_all=False)
print(len(fofmass))
print(fofpos.shape)
maxpos=numpy.where(fofmass==max(fofmass))
print(maxpos)
fofpos=fofpos[maxpos,:]
print(fofpos)
'''
'''
fig, ax = plt.subplots(1,1,figsize=(11,9))
mets = []
reds = []
ids=[]
fofmass=[]
fofpos=[]
zs   = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
for red in zs:
    #metallicity,output_redshift=line_plots(basePath,red,0,'Metallicity')
    #bh_IDs,output_redshift=arepo_package.get_particle_property(basePath,'ParticleIDs',5,red,list_all=True)
    fofmassin,output_redshift=arepo_package.get_group_property(basePath,'GroupMass',red,list_all=True)
    fofposin,output_redshift=arepo_package.get_group_property(basePath,'GroupPos',red,list_all=False)
    #metallicity=gaussian_filter1d(metallicity,sigma=1)
    print(output_redshift)
    #mets.append(numpy.nanmean(metallicity))
    #ids.append(len(bh_IDs))
    fofmass.append(max(fofmassin)*1e10)
    fofpos.append(fofposin)
    reds.append(output_redshift)

maxpos=numpy.where(fofmass==max(fofmass))
#print reds
#print fofmass
#plt.plot(reds,numpy.log10(mets))
plt.plot(reds,numpy.log10(fofmass))
plt.axhline(y=numpy.log10(5e10),color='r')
#ax.set_ylabel('log gas metallicity $(Fe/H)/(Fe/H)_{\odot}$',fontsize=15)
ax.set_ylabel('log FOF Mass ($M_{\odot}/h$)',fontsize=15)
ax.set_xlabel('redshift',fontsize=15)
plt.title('FOF Mass Level Max 9',size=15)
fig.savefig('Prof_Package_Test')
'''
'''
#-----Warning: Sublink merger trees have been computed only up to z=5. DO NOT select a root redshift less than 5, for now------------------------------------- 
root_subhalo_index=1   
root_redshift=19
#------------------------Get the indices and snapshots of the most massive progenitors---------------------------
Progenitor_SubhaloIndices,Progenitor_Snaps=arepo_package.get_sublink_progenitors_most_massive_branch(basePath,root_subhalo_index,root_redshift)
#-----------------------------------------------------------------------------------------------------------------

#------------------------Select the oldest progenitor and trace its descendants. -----------------------------------
#------------------------Check if the descendant track matches the most massive progenitor track--------------------
oldest_progenitor_subhalo_index=Progenitor_SubhaloIndices[-1]
oldest_progenitor_snap=Progenitor_Snaps[-1]
snap_list,redshift_list=arepo_package.get_snapshot_redshift_correspondence(basePath)
oldest_progenitor_redshift=redshift_list[snap_list==oldest_progenitor_snap][0]
Descendant_SubhaloIndices,Descendant_Snaps=arepo_package.get_sublink_descendants(basePath,oldest_progenitor_subhalo_index,oldest_progenitor_redshift)
Descendant_SubhaloIndices_upto_progenitor=Descendant_SubhaloIndices[0:len(Progenitor_SubhaloIndices)]

print("Progenitor indices are:",Progenitor_SubhaloIndices)
print("All Descendant indices:",Descendant_SubhaloIndices)
print("All Descendant indices upto the progenitor snapshot in reverse:",Descendant_SubhaloIndices_upto_progenitor[::-1])
print ("Difference (If they are all 0, it is good. If not, there's a problem, report to Aklant):",Progenitor_SubhaloIndices-Descendant_SubhaloIndices_upto_progenitor[::-1])
'''

