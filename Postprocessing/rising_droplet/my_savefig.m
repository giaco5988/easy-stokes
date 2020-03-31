%save figures in a folder in a certian file extension

close all
clear all
path = '~/Documents/Phd_projects/Viscous_drop_buoyancy/Report9_flow_field/figures/';

%phisical variable
lambda = [0.5 5];
Ca = 6;

for i=1:numel(lambda)
    for k = 1:numel(Ca)
        
        if lambda(i) == 0.5
            
            name = 'lambda05_Ca6_opt_oblate1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda05_Ca6_opt_oblate_vorticity1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda05_Ca6_opt_oblate2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda05_Ca6_opt_oblate_vorticity2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda05_Ca6_opt_prolate1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda05_Ca6_opt_prolate_vorticity1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda05_Ca6_opt_prolate2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda05_Ca6_opt_prolate_vorticity2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
        elseif lambda(i) == 5
            
            name = 'lambda5_Ca6_opt_oblate1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda5_Ca6_opt_oblate_vorticity1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda5_Ca6_opt_oblate2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda5_Ca6_opt_oblate_vorticity2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda5_Ca6_opt_prolate1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda5_Ca6_opt_prolate_vorticity1';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda5_Ca6_opt_prolate2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
            name = 'lambda5_Ca6_opt_prolate_vorticity2';
            open([name '.fig'])
            print('-depsc','-r300',[path name '.eps'])
            
        end
        
    end
end

close all

%fix_lines(['./Figures/h_t_theta_A.eps'],['./Figures/h_t_theta_A.eps'])