%compute interface shape for droplet in extensional flow usinf Newton
%method

function tutorial_newtonMethodSpectralXYmodesFunction(PARAM)

%current directory
here = pwd;

%print to screen
PARAM.ODE = 0;
printToScreen(PARAM);

%initial condition
[initial,x,y] = initialConditionDrop(PARAM);

%coordinates for plotting
aIN = x;  bIN = y;

%dealiasing
xBase = initial(1:2:end-1);
yBase = initial(2:2:end);
xBase = xBase(1:PARAM.dealiasing);
yBase = yBase(1:PARAM.dealiasing);

%final volume that I want
V0 = PARAM.V0;

%position of the center of mass (needed only if in drop frame)
xcm = CenterMassCurvAxisSpectral(x,y,PARAM);

%initialize
manyRES = zeros(1,PARAM.stop);
if PARAM.dropFrame==0
    perturb = zeros(PARAM.dealiasing-1,1);
elseif PARAM.dropFrame==1
    perturb = zeros(PARAM.dealiasing-2,1);
end

%nonlinear equation
fNonlinear = @(unk) tutorial_NormalVelocitySpectralVolumeXY2modes(unk,xBase,yBase,V0,xcm,PARAM);
[~,~,ny] = fNonlinear(perturb);

%check initial shape
if min(ny)<-1e-5
    warning('Normal vector of the initial shape is not accurate')
end

% Newton iterations
disp('Newton iteration starts')
    
quit=0;count=0;
while ~quit
    
    %NON LINEAR SOLUTION
    %normal velocity
    [u,~,~,x,y,firstMode,K1,K2] = fNonlinear(perturb);
    R = u;
    
    %Jacobian
    J = JacobianHandle(fNonlinear,perturb,PARAM.dh);
    
    %compute ellipticity
    L = max(x)-min(x);  B = 2*y(round(numel(y)/2));
    D = (L-B)/(L+B);
    
    %new guess
    step = J\R;
    perturb = perturb - step*PARAM.CutStep;
    
    % convergence test
    res = max(abs([R; firstMode']));
    %res = norm([R; firstMode]);
    manyRES(count+1) = res;
    figure(1)
    plot([x; flip(x)],[y; -flip(y)],[aIN; flip(aIN)],[bIN; -flip(bIN)],'--');
    grid on; axis equal; axis([-4 4 -2 2]); drawnow;
    disp(['Iteration ' num2str(count+1) ' res=' num2str(res) '  D=' num2str(D)]);
    legend('Solution','Initial guess')
    title(['Droplet shape Ca' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc)])
    if count>PARAM.stop||res>1e1||isnan(res)==1
        disp('Limit number of iteration reached, stop!')
        if res>1
            disp('no convergence');
        end
        break;
    end
    if res < PARAM.ResBreak; quit=1; disp('converged'); continue; end
    
    %count newton iteration
    count = count+1;
    
end

if PARAM.plotRes==1
    manyRES = manyRES(1:count+1);
    figure(4)
    semilogy(1:count+1,manyRES,'-o')
    hold on
    grid on
    xlabel('$n$','interpreter','latex')
    ylabel('$||\hat{\mathbf{u}}^{(n)}||_\infty$','interpreter','latex')
    title('Convergence')
end

if PARAM.plotCurv==1
    %plot curvature
    figure(5)
    plot(PARAM.t,K1+K2,'-')
    grid on
    hold on
    %plot(PARAM.t,K2,'-')
    %plot(PARAM.t,K1+K2,'-k')
    title(['Curvature Ca=' num2str(PARAM.Ca)])
    %legend('K_1','K_2','K_1+K_2','Location','Best')
    xlabel('$s$','interpreter','latex')
    ylabel('$K$','interpreter','latex')
end

cd(PARAM.res)
save(PARAM.filename)
cd(here)