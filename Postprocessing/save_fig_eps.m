%save current figure in .eps to a specified folder and come back to the current
%folder

Gt = 0.1;
figname = ['aft',num2str(100*Gt)];
open(strcat(figname,'.fig'))

filename = '~/Documents/Phd_projects/Viscous_drop_buoyancy/Report5_stab_unstab/figures';
cd(filename)

filename = strcat(figname,'.eps');
%problem: the exported files are in black and white!!!
saveas(gcf,filename)
%export_fig ciao.eps
%print(gcf,filename)

cd ~/Documents/MATLAB/droplet_simulations/results/delta_graph/

close all