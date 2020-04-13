import sys
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
#sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
sys.path.append('/home/kkosciw/Python_Scripts/arepo_package')
import site
site.addsitedir("/home/kkosciw/Python_Scripts")


import arepo_package
#%pylab inline

#from arepo_script import *

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import astro_constants as ac


def profile_plots_2D(basePath_uniform,desired_particle_property,p_type,desired_redshift,index_of_selected_halo,projection_plane,slice_thickness,bin_number,title_name,save_name):
    #particle_property can equal Density, Temperature, Metallicity
    #projection_plane will be xy,xz,or yz
    
    boxsize=arepo_package.get_box_size(basePath_uniform)
    desired_redshift_of_selected_halo=desired_redshift
 
    particle_property='Coordinates'
    p_type=0
    gas_particle_positions,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo,group_type='groups',list_all=True)

    xpos=gas_particle_positions[:,0]
    ypos=gas_particle_positions[:,1]
    zpos=gas_particle_positions[:,2]

    #First Mask
    if projection_plane=='xy':
        left=numpy.average(zpos)-(slice_thickness/2.0)
        right=numpy.average(zpos)+(slice_thickness/2.0)
        mask1=(zpos>left)&(zpos<right)
    if projection_plane=='xz':
        left=numpy.average(ypos)-(slice_thickness/2.0)
        right=numpy.average(ypos)+(slice_thickness/2.0)
        mask1=(ypos>left)&(ypos<right)
    if projection_plane=='yz':
        left=numpy.average(xpos)-(slice_thickness/2.0)
        right=numpy.average(xpos)+(slice_thickness/2.0)
        mask1=(xpos>left)&(xpos<right)
 
    xpos=xpos[mask1]
    ypos=ypos[mask1]
    zpos=zpos[mask1]

    if (desired_particle_property=='Metallicity'):
        particle_property='GFM_Metallicity'
        gas_metallicity,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo,group_type='groups',list_all=False)
        gas_metals=gas_metallicity[mask1]/0.0127
        particle_prop=gas_metals

    if (desired_particle_property=='Temperature'):
        particle_property='InternalEnergy'
        particle_internal_energy,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo,group_type='groups',list_all=False)
        particle_property='ElectronAbundance'
        particle_electron_abundance,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo,group_type='groups',list_all=False)
        g_minus_1 = (5.0/3.0) - 1.0
        XH = 0.76
        mu=(4*ac.MP)/(1+3*XH+4*XH*particle_electron_abundance[mask1])
        gas_temperature = g_minus_1*(particle_internal_energy[mask1]/ac.KB)*(10**10)*mu
        particle_prop=gas_temperature

    if (desired_particle_property=='Density'):
        particle_property='Masses'
        particle_mass,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo,group_type='groups',list_all=False) 
        particle_prop=particle_mass[mask1]*1e10
   
    ndx=bin_number
    ndy=bin_number
    ndz=bin_number

    x=numpy.linspace(min(xpos),max(xpos),ndx)
    y=numpy.linspace(min(ypos),max(ypos),ndy)
    z=numpy.linspace(min(zpos),max(zpos),ndz)
    dx=numpy.diff(x)[0]
    dy=numpy.diff(y)[0]
    dz=numpy.diff(z)[0]
    proj_property=[]

    #2nd Mask
    if (projection_plane=='xy'):
        First,Second=numpy.meshgrid(x,y)
        for yi in y:
            for xi in x:
                mask_x=(xpos>(xi-dx/2.0))&(xpos<(xi+dx/2.))
                mask_y=(ypos>(yi-dy/2.))&(ypos<(yi+dy/2.))
                mask2=(mask_x) & (mask_y)
                if (desired_particle_property=='Density'):
                    proj_property.append(numpy.sum(particle_prop[mask2])/(dx*dy*dz))
                else:
                    proj_property.append(numpy.average(particle_prop[mask2]))
    if (projection_plane=='xz'):
        First,Second=numpy.meshgrid(x,z)
        for zi in z:
            for xi in x:
                mask_x=(xpos>(xi-dx/2.0))&(xpos<(xi+dx/2.))
                mask_z=(zpos>(zi-dz/2.))&(zpos<(zi+dz/2.))
                mask2=(mask_x) & (mask_z)
                if (desired_particle_property=='Density'):
                    proj_property.append(numpy.sum(particle_prop[mask2])/(dx*dy*dz))
                else:
                    proj_property.append(numpy.average(particle_prop[mask2]))
    if (projection_plane=='yz'):
        First,Second=numpy.meshgrid(y,z)
        for zi in z:
            for yi in y:
                mask_z=(zpos>(zi-dz/2.0))&(zpos<(zi+dz/2.))
                mask_y=(ypos>(yi-dy/2.))&(ypos<(yi+dy/2.))
                mask2=(mask_z) & (mask_y)
                if (desired_particle_property=='Density'):
                    proj_property.append(numpy.sum(particle_prop[mask2])/(dx*dy*dz))
                else:
                    proj_property.append(numpy.average(particle_prop[mask2]))
    
    proj_property=numpy.asarray(proj_property)
    proj_property[numpy.isnan(proj_property)]=1e-19

    Proj_property=proj_property.reshape(bin_number,bin_number)

    fig = plt.figure()
    ax = fig.add_subplot(111)
   

    if (desired_particle_property=='Density'):
        plt.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(particle_prop)/(dx*dy*dz),vmax=Proj_property.max()))
        cbar=plt.colorbar()
        cbar.set_label('log gas density ($M_{\odot}/(cKpc/h)^3$)')
    if (desired_particle_property=='Metallicity'):
        plt.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(particle_prop),vmax=Proj_property.max()))
        cbar=plt.colorbar()
        cbar.set_label('log gas metallicity $(Fe/H)/(Fe/H)_{\odot})$')
    if (desired_particle_property=='Temperature'):
        plt.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(particle_prop),vmax=Proj_property.max()))
        cbar=plt.colorbar()
        cbar.set_label('log gas temperature $K$')

    plt.title(title_name)
    fig.savefig(save_name)
    return 0



path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'

uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
#uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass5.90_logFOFseedmass10.70/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'

profile_plots_2D(basePath_uniform,'Density',0,0.2,0,'xy',100,100,'Density z=0.2','Prof_Package_Test')
