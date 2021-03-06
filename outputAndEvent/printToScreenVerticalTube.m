%print to screen information when starting the simulation

function printToScreenVerticalTube(L,x0,alpha,Bond,lambda,dropFrame,Ca)

%output hydrodynamics
disp('HYDRODYNAMIC PARAMETERS')
disp(['Bo=' num2str(Bond) ', Ca=' num2str(Ca) ', lambda=' num2str(lambda)])

%ouput geometry
disp('GEOMETRY')
disp(['Channel lenght is L=' num2str(L)])

%ouput bubble
disp('BUBBLE GEOMETRY')
disp(['Initial bubble position is x0=' num2str(x0) ' with unperturbed radius r0=' num2str(alpha)])

%frame of reference
if dropFrame==0
    frame = 'lab frame';
elseif dropFrame==1
    frame = 'drop frame';
elseif dropFrame==2
    frame = 'drop frame but replace drop in center of mass';
else
    error('Not implemented')
end
disp(['Simulation is in the ' frame])

disp(['The capillary lenght is lc=' num2str(sqrt(1/Bond)) ', when lc>1.1 (rouglhy) the droplet might get stucked'])