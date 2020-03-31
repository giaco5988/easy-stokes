%open picture .fig and save it in .eps

clear variables
close all

source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/figures';
dest = '~/Documents/phd_projects/micromotor/report5_no_injection/figures/';
cd(source);

alpha = 1.5:0.1:1.8;

for i = 1:numel(alpha)
    
    disp(i)
    
    %open figure
    name = ['alpha' num2str(alpha(i)*10)];
    open([name '.fig'])
    
    %save in eps to folder
    cd(dest)
    print(name,'-depsc')
    cd(source)
     
end