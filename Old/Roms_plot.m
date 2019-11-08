inputdata = 'ROMS_500m_ur_and_ud_VSF_W_binned.mat';
datapath ='/Volumes/bumpusdata/Impact-of-Convergence-Zones-on-Lagrangian-Structure-Function-Statistics-in-the-Gulf-of-Mexico/Data/Eulerian/';


load([datapath inputdata])

loglog(r,dsf,r,rsf)
legend('div','rot')