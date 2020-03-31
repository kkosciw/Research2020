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

import astro_constants as ac
'''
#Gas Density

min_edge=-1
max_edge=3
Nbins=50
desired_redshift_of_selected_halo=4

p_type=0  #0 is for gas particles
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1005831581
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_IDs)
all_bh_masses,output_redshift=arepo_package.get_particle_property(basePath_uniform,'Masses',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_masses)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]

FOF_bh_IDs,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,group_type='groups',list_all=False)
print(FOF_bh_IDs)
FOF_bh_masses,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,'Masses',5,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,group_type='groups',list_all=False)
print(FOF_bh_masses)

bin_centers10,mass_distribution10,mass_density10=arepo_package.get_halo_density_profile(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)



p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_IDs)

all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)

index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]

FOF_bh_IDs,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,group_type='groups',list_all=False)
print(FOF_bh_IDs)
FOF_bh_masses,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,'Masses',5,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,group_type='groups',list_all=False)
print(FOF_bh_masses)

bin_centers9,mass_distribution9,mass_density9=arepo_package.get_halo_density_profile(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)


p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
#seed_bh_mass=8e5
#min_FOF_seed_mass=5e10
uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass5.90_logFOFseedmass10.70/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1056102916
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_IDs)
all_bh_masses,output_redshift=arepo_package.get_particle_property(basePath_uniform,'Masses',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_masses)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]

FOF_bh_IDs,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,group_type='groups',list_all=False)
print(FOF_bh_IDs)
FOF_bh_masses,output_redshift=arepo_package.get_particle_property_within_groups(basePath_uniform,'Masses',5,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,group_type='groups',list_all=False)
print(FOF_bh_masses)

bin_centers11,mass_distribution11,mass_density11=arepo_package.get_halo_density_profile(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)




#f,ax=plt.subplots(1,1,figsize=(10,8))
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_centers9,mass_density9*1e10,label='levelmin=7,levelmax=9',color='red')
plt.plot(bin_centers10,mass_density10*1e10,label='levelmin=7,levelmax=10',color='blue')
plt.plot(bin_centers11,mass_density11*1e10,label='levelmin=7,levelmax=11',color='green')
plt.yscale('log')
#ax.tick_params(labelsize=30)
plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='lower left')#,fontsize=25)

plt.ylabel('gas density ($M_{\odot}/(cKpc/h)^3$)')#,fontsize=30)

#plt.ylim(1e-1,1e8)
plt.title('Gas Density at z=1')
fig.savefig('Gas_Density_1D_1_plot')
'''
'''
#Gas Metallicity

min_edge=-1
max_edge=3
Nbins=50

#Level Max 10
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
desired_redshift_of_selected_halo=0
p_id=1005831581
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_IDs)
all_bh_masses,output_redshift=arepo_package.get_particle_property(basePath_uniform,'Masses',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_masses)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]

bin_centers10,mass_distribution10,mass_density10,metallicity_distribution10,velocity_distribution10,temperature10,distance_distribution10=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)



#Level Max 9
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_IDs)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers9,mass_distribution9,mass_density9,metallicity_distribution9,velocity_distribution9,temperature9,distance_distribution9=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)

#Level Max 11
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
seed_bh_mass=8e5
min_FOF_seed_mass=5e10
uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass%.2f_logFOFseedmass%.2f/AREPO'%(numpy.log10(seed_bh_mass),numpy.log10(min_FOF_seed_mass))
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1056102916
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_IDs)
all_bh_masses,output_redshift=arepo_package.get_particle_property(basePath_uniform,'Masses',5,desired_redshift_of_selected_halo,list_all=False)
print(all_bh_masses)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers11,mass_distribution11,mass_density11,metallicity_distribution11,velocity_distribution11,temperature11,distance_distribution11=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)


#Plot the profile
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_centers9,metallicity_distribution9/0.0127,label='levelmin=7,levelmax=9',color='red')
plt.plot(bin_centers10,metallicity_distribution10/0.0127,label='levelmin=7,levelmax=10',color='blue')
plt.plot(bin_centers11,metallicity_distribution11/0.0127,label='levelmin=7,levelmax=11',color='green')

plt.yscale('log')

plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='lower left')#,fontsize=25)

plt.ylabel('gas metallicity $(Fe/H)/(Fe/H)_{\odot})$')#,fontsize=30)
plt.title('Gas Metallicity at z=3')
fig.savefig('Metallicity_1D_3_plot')
'''


'''
#Gas Velocity Magnitude

min_edge=-1
max_edge=3
Nbins=50

#Level Max 10
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
desired_redshift_of_selected_halo=3
p_id=1005831581
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]

bin_centers10,mass_distribution10,mass_density10,metallicity_distribution10,velocity_distribution10,temperature10,distance_distribution10=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)

#Level Max 9
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers9,mass_distribution9,mass_density9,metallicity_distribution9,velocity_distribution9,temperature9,distance_distribution9=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)


#Level Max 11
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
seed_bh_mass=8e5
min_FOF_seed_mass=5e10
uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass%.2f_logFOFseedmass%.2f/AREPO'%(numpy.log10(seed_bh_mass),numpy.log10(min_FOF_seed_mass))
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1056102916
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers11,mass_distribution11,mass_density11,metallicity_distribution11,velocity_distribution11,temperature11,distance_distribution11=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)



#Plot the profile
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_centers9,velocity_distribution9,label='levelmin=7,levelmax=9',color='red')
plt.plot(bin_centers10,velocity_distribution10,label='levelmin=7,levelmax=10',color='blue')
plt.plot(bin_centers11,velocity_distribution11,label='levelmin=7,levelmax=11',color='green')

plt.yscale('log')

plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='upper left')#,fontsize=25)

plt.ylabel('gas velocity magnitude $km \sqrt{a}/s$')#,fontsize=30)
plt.title('Gas Velocity z=3')
fig.savefig('Velocity_1D_3_plot')
'''
'''
#Dark Matter Density


min_edge=-1
max_edge=3
Nbins=50


#Level Max 10
p_type=1
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
desired_redshift_of_selected_halo=4
#path_to_output='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
#run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
#basePath=path_to_output+run+'/output/'
print('Header file contents')
header=arepo_package.load_snapshot_header(basePath_uniform, 0.5)
MassTable=header['MassTable']
print("MassTable for all particles:",MassTable)
particle_type1_mass=MassTable[1]
print("Mass of high resolution DM particle in M_sun:",particle_type1_mass*(10**10))

p_id=1005831581
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]

bin_centers10,distance_distribution10=arepo_package.get_DM_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)

DM_density10 = distance_distribution10*particle_type1_mass*(10**10)*3/4/3.14/(10**bin_centers10)**3

#print(DM_density10)
#print(distance_distribution10)
#print(len(distance_distribution10))
#print(bin_centers10)
#print(len(bin_centers10))
#print(particle_type1_mass)


#Level Max 9
p_type=1
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983
print('Header file contents')
header=arepo_package.load_snapshot_header(basePath_uniform, 0.5)
MassTable=header['MassTable']
print("MassTable for all particles:",MassTable)
particle_type1_mass=MassTable[1]
print("Mass of high resolution DM particle in M_sun:",particle_type1_mass*(10**10))
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers9,distance_distribution9=arepo_package.get_DM_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)

DM_density9 = distance_distribution9*particle_type1_mass*(10**10)*3/4/3.14/(10**bin_centers9)**3


#Level Max 11
p_type=1
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
seed_bh_mass=8e5
min_FOF_seed_mass=5e10
uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass%.2f_logFOFseedmass%.2f/AREPO'%(numpy.log10(seed_bh_mass),numpy.log10(min_FOF_seed_mass))
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1056102916
print('Header file contents')
header=arepo_package.load_snapshot_header(basePath_uniform, 0.5)
MassTable=header['MassTable']
print("MassTable for all particles:",MassTable)
particle_type1_mass=MassTable[1]
print("Mass of high resolution DM particle in M_sun:",particle_type1_mass*(10**10))
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers11,distance_distribution11=arepo_package.get_DM_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)

DM_density11 = distance_distribution11*particle_type1_mass*(10**10)*3/4/3.14/(10**bin_centers11)**3



#Plot the profile
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_centers9,DM_density9,label='levelmin=7,levelmax=9',color='red')
plt.plot(bin_centers10,DM_density10,label='levelmin=7,levelmax=10',color='blue')
plt.plot(bin_centers11,DM_density11,label='levelmin=7,levelmax=11',color='green')

plt.yscale('log')

plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='lower left')#,fontsize=25)

plt.ylabel('DM Density ($M_{\odot}/h)/(cKpc/h)^3$')#,fontsize=30)
plt.title('DM Density z=4')
fig.savefig('DM_Density_1D_4_plot')
'''
'''
#Gas Temperature

min_edge=-1
max_edge=3
Nbins=50

#Level Max 10
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax10_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
desired_redshift_of_selected_halo=7
p_id=1005831581
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)

all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]

bin_centers10,mass_distribution10,mass_density10,metallicity_distribution10,velocity_distribution10,temperature10,distance_distribution10=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)


#Level Max 9
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers9,mass_distribution9,mass_density9,metallicity_distribution9,velocity_distribution9,temperature9,distance_distribution9=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)

#Level Max 11
p_type=0
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
seed_bh_mass=8e5
min_FOF_seed_mass=5e10
uniform_run='L25n128MUSIC_rerun_zoom_levelmax11_haloindex100_redshift0.00_logbhseedmass%.2f_logFOFseedmass%.2f/AREPO'%(numpy.log10(seed_bh_mass),numpy.log10(min_FOF_seed_mass))
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1056102916
all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,desired_redshift_of_selected_halo,list_all=False)
all_group_ids=arepo_package.generate_group_ids(basePath_uniform,desired_redshift_of_selected_halo,5,save_output_path='./',create=True)
index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
bin_centers11,mass_distribution11,mass_density11,metallicity_distribution11,velocity_distribution11,temperature11,distance_distribution11=arepo_package.get_halo_profile_property(basePath_uniform,p_type,desired_redshift_of_selected_halo,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)



#Plot the profile
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_centers9,temperature9,label='levelmin=7,levelmax=9',color='red')
plt.plot(bin_centers10,temperature10,label='levelmin=7,levelmax=10',color='blue')
plt.plot(bin_centers11,temperature11,label='levelmin=7,levelmax=11',color='green')

plt.yscale('log')

plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='upper right')#,fontsize=25)

plt.ylabel('gas temperature $K$')#,fontsize=30)
plt.title('Temperature at z=7')
fig.savefig('Temperature_1D_7_plot')
'''
'''

#Gas Density Multiple Red Shifts
min_edge=-1
max_edge=3
Nbins=50
desired_redshift_of_selected_halo=[0,0.5,1,2,3]

p_type=0  #0 is for gas particles
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983

bin_center=[]
mass_densities=[]

for i in desired_redshift_of_selected_halo:
    print(len(bin_center))
    all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,i,list_all=False)
    all_bh_masses,output_redshift=arepo_package.get_particle_property(basePath_uniform,'Masses',5,i,list_all=False)
    all_group_ids=arepo_package.generate_group_ids(basePath_uniform,i,5,save_output_path='./',create=True)
    index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
    bin_centers,mass_distribution,mass_density=arepo_package.get_halo_density_profile(basePath_uniform,p_type,i,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)
    bin_center=numpy.append(bin_center,bin_centers)
    mass_densities=numpy.append(mass_densities,mass_density)
    


fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_center[0:48],mass_densities[0:48]*1e10,label='z=0',color='red')
plt.plot(bin_center[49:98],mass_densities[49:98]*1e10,label='z=0.5',color='green')
plt.plot(bin_center[99:148],mass_densities[99:148]*1e10,label='z=1',color='blue')
plt.plot(bin_center[149:198],mass_densities[149:198]*1e10,label='z=2',color='magenta')
plt.plot(bin_center[199:248],mass_densities[199:248]*1e10,label='z=3',color='orange')
#plt.plot(bin_center[250:299],mass_densities[250:299]*1e10,label='z=4',color='cyan')

plt.yscale('log')

plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='lower left')#,fontsize=25)

plt.ylabel('gas density, ($M_{\odot}/(cKpc/h)^3$)')#,fontsize=30)

plt.title('Gas Density Level Max 9')
fig.savefig('Gas_Density_Level_Max_9')
'''
'''
#Velocity,Temperature,Metallicity,Stellar Density Multiple Redshifts

min_edge=-1
max_edge=3
Nbins=50
desired_redshift_of_selected_halo=[0,0.5,1,2,3]

p_type=4  #0 is for gas particles, 4 is for star particles
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983

bin_center=[]
prop=[]

for i in desired_redshift_of_selected_halo:
    print(len(bin_center))
    all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,i,list_all=False)
    all_group_ids=arepo_package.generate_group_ids(basePath_uniform,i,5,save_output_path='./',create=True)
    index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
    bin_centers,mass_distribution,mass_density,metallicity_distribution,velocity_distribution,distance_distribution=arepo_package.get_halo_profile_property(basePath_uniform,p_type,i,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)
    #Add temperature to previous line when done with stellar density
    bin_center=numpy.append(bin_center,bin_centers)
    prop=numpy.append(prop,mass_distribution*1e10*3/4/3.14/(10**bin_centers)**3)



fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_center[0:48],prop[0:48],label='z=0',color='red')
plt.plot(bin_center[49:98],prop[49:98],label='z=0.5',color='green')
plt.plot(bin_center[99:148],prop[99:148],label='z=1',color='blue')
plt.plot(bin_center[149:198],prop[149:198],label='z=2',color='magenta')
plt.plot(bin_center[199:248],prop[199:248],label='z=3',color='orange')
#plt.plot(bin_center[249:298],prop[249:298]*1e10,label='z=4',color='cyan')

plt.yscale('log')
plt.ylim(1e-1,1e10)

plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='lower left')#,fontsize=25)

#plt.ylabel('gas temperature $K$')#,fontsize=30)
#plt.ylabel('gas velocity magnitude $km \sqrt{a}/s$')#,fontsize=30)
#plt.ylabel('gas metallicity $(Fe/H)/(Fe/H)_{\odot})$')#,fontsize=30)
plt.ylabel('Stellar Density ($M_{\odot}/h)/(cKpc/h)^3$')#,fontsize=30)



plt.title('Stellar Density Level Max 9')
fig.savefig('Stellar_Density_Level_Max_9')
'''
'''
#DM Density Multiple Redshifts


min_edge=-1
max_edge=3
Nbins=50
desired_redshift_of_selected_halo=[0,0.5,1,2,3]

p_type=1
path_to_uniform_run='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
uniform_run='L25n128MUSIC_rerun_zoom_levelmax9_haloindex100_redshift0.00/AREPO'
basePath_uniform=path_to_uniform_run+uniform_run+'/output_BH_NGB_256/'
p_id=1002600983


print('Header file contents')
header=arepo_package.load_snapshot_header(basePath_uniform, 0.5)
MassTable=header['MassTable']
print("MassTable for all particles:",MassTable)
particle_type1_mass=MassTable[1]
print("Mass of high resolution DM particle in M_sun:",particle_type1_mass*(10**10))

bin_center=[]
dm_densities=[]

for i in desired_redshift_of_selected_halo:
    all_bh_IDs,output_redshift=arepo_package.get_particle_property(basePath_uniform,'ParticleIDs',5,i,list_all=False)
    all_group_ids=arepo_package.generate_group_ids(basePath_uniform,i,5,save_output_path='./',create=True)
    index_of_selected_halo_zoom=(all_group_ids[all_bh_IDs==p_id])[0]
    bin_centers,distance_distribution=arepo_package.get_DM_property(basePath_uniform,p_type,i,index_of_selected_halo_zoom,min_edge,max_edge,Nbins,CENTER_AROUND='MOST_MASSIVE_BLACKHOLE',p_id=p_id)
    DM_density =distance_distribution* particle_type1_mass*(10**10)*3/4/3.14/(10**bin_centers)**3
    bin_center=numpy.append(bin_center,bin_centers)
    dm_densities=numpy.append(dm_densities,DM_density)


fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(bin_center[0:48],dm_densities[0:48],label='z=0',color='red')
plt.plot(bin_center[49:98],dm_densities[49:98],label='z=0.5',color='green')
plt.plot(bin_center[99:148],dm_densities[99:148],label='z=1',color='blue')
plt.plot(bin_center[149:198],dm_densities[149:198],label='z=2',color='magenta')
plt.plot(bin_center[199:248],dm_densities[199:248],label='z=3',color='orange')
#plt.plot(bin_center[249:298],prop[249:298]*1e10,label='z=4',color='cyan')

plt.yscale('log')

plt.xlabel('log10(distance from halo center) in kpc/h')#,fontsize=30)

plt.legend(loc='lower left')#,fontsize=25)

plt.ylabel('DM Density ($M_{\odot}/h)/(cKpc/h)^3$')#,fontsize=30)

plt.title('DM Density Level Max 9')
fig.savefig('DM_Density_Level_Max_9')
'''


