#!/usr/bin/env python
# coding: utf-8

# # Mesogens with NP | Equilibrium

# ## Temperature 7.2 | GPU Run

# ### Date: 12/10/2019 | System P = 1.8, Expected value of $T_c$ : 7.09 |

# In[ ]:


from __future__ import division
import hoomd
import hoomd.md


# In[ ]:


#-----Define relevant variables
p_max = 1.8;
t_max = 7.2;
copies = 1;
steps_run = 1e5;
init_file = "T_CM&NP_" + str(t_max) + "_P_" + str(p_max) + "_ramp.gsd"


# In[ ]:


#-----Define a simulation context

    #-----This is if I want to run on a GPU

#hoomd.context.initialize("--mode=gpu");

    #-----This is if I want to run on a CPU

hoomd.context.initialize("--mode=cpu");

# In[ ]:


#-----Extract the configuration of the system and expand the system

snap = hoomd.data.gsd_snapshot(init_file, frame = -1);
snap.replicate(copies,copies,copies);
system = hoomd.init.read_snapshot(snap);


# In[ ]:


#-----Define each mesogen in the local reference frame of each center of mass
rigid = hoomd.md.constrain.rigid();
rigid.set_param('M', 
               types = ['A']*8,
               positions = [(-4,0,0),(-3,0,0),(-2,0,0),(-1,0,0),
                            (1,0,0),(2,0,0),(3,0,0),(4,0,0)]);


# In[ ]:


#-----Declare molecules as rigid bodies
rigid.create_bodies();


# In[ ]:


#-----Define the potential energy
nl = hoomd.md.nlist.tree();
lj = hoomd.md.pair.lj(r_cut = 3.5, nlist = nl);
lj.set_params(mode = 'shift')


# In[ ]:


#------Define the interaction
lj.pair_coeff.set('NP','NP', epsilon = 1.0, sigma = 5.0);
lj.pair_coeff.set('M', 'M', epsilon = 1.0, sigma = 1.0);
lj.pair_coeff.set('A', 'A', epsilon = 1.0, sigma = 1.0);
lj.pair_coeff.set('M', 'A', epsilon = 1.0, sigma = 1.0);
lj.pair_coeff.set('NP', 'M', epsilon = 1.0, sigma = 3.0);
lj.pair_coeff.set('NP', 'A', epsilon = 1.0, sigma = 3.0);


# In[ ]:


#------Select an standar integrator
hoomd.md.integrate.mode_standard(dt = 0.005);

#-----Define some groups and make their union

nanoparticles = hoomd.group.type(name = 'NPs', type = 'NP');
mesogens = hoomd.group.rigid_center();
groupNP_mes = hoomd.group.union(name = 'NP_Mes', a = nanoparticles, b = mesogens);


# In[ ]:


#-----Integrate using NPT

npt = hoomd.md. integrate.npt(group = groupNP_mes, kT = t_max, tau = 14.1, tauP = 14.1, P = p_max);


# In[ ]:


#-----Save data

log_file = "T_" + str(t_max) + "_P_" + str(p_max) + "_equilibrium.log"
gsd_file = "T_" + str(t_max) + "_P_" + str(p_max) + "_equilibrium.gsd"
meso_gsd_file = "T_CM_" + str(t_max) + "_P_" + str(p_max) + "_equilibrium.log"

log = hoomd.analyze.log(filename = log_file,
                       quantities = ['num_particles', 
                                    'ndof',
                                    'translational_ndof',
                                    'rotational_ndof',
                                    'potential_energy',
                                    'kinetic_energy',
                                    'translational_kinetic_energy',
                                    'rotational_kinetic_energy',
                                    'temperature',
                                    'pressure',
                                    'volume'],
                       period = 1e2,
                       overwrite = True);
gsd = hoomd.dump.gsd(gsd_file, period = 1e2, group = hoomd.group.all(), overwrite = True);
meso_gsd = hoomd.dump.gsd(meso_gsd_file, period=1e2, group = mesogens, overwrite = True);


# In[ ]:


#-----Run the simulation

hoomd.run(steps_run)


# In[ ]:


#-----Get volume and density information.
system.box.get_volume()


# In[ ]:


system.get_metadata()

