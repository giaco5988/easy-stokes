function stop = outfun_newton(x,optimValues,state,nonLinearFunction)

%% COMPUTE NONLINEAR FUNCTION (GETS CURRENT SHAPE)
[~,~,~,xGrid,yGrid,~,xBase,yBase] = nonLinearFunction(x);

%% PLOT CURRENT SHAPE
if state == 'iter'
    figure(1)
    plot([xGrid, flip(xGrid)],[yGrid, -flip(yGrid)])
    hold on
    axis equal
    plot([xBase; flip(xBase)],[yBase; -flip(yBase)],'--')
    legend('Current value','Initial guess')
    axis([-3 3 -3 3])
    xlabel('z')
    ylabel('r')
    hold off
end

%% PLOT RESIDUALS
if state == 'iter'
    figure(2)
    res = norm(optimValues.fval,inf);
    semilogy(optimValues.iteration,res,'or','MarkerSize',10)
    hold on
    semilogy(optimValues.iteration,res,'.k','MarkerSize',30)
    xyLabelTex('\mathrm{iteration}','||R||_\infty')
    grid on
end

%% PLOT CURVATURE OF THE FINAL SHAPE
if state=='done'
    figure(3)
    [k1,k2] = computeCurvatureSplines(xGrid,yGrid,1);
    k = k1+k2;
    plot(xGrid,k)
    xyLabelTex('z','K')
    title('Curvature')
    grid on
end

stop = 0;