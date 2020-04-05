function stop = outfun_newton(x,optimValues,state,nonLinearFunction)

%% COMPUTE NONLINEAR FUNCTION (GETS CURRENT SHAPE)
[~,~,~,xBaseGrid,yBaseGrid,~,xBase,yBase] = nonLinearFunction(x);

%% PLOT CURRENT SHAPE
if state == 'iter'
    figure(1)
    plot([xBaseGrid, flip(xBaseGrid)],[yBaseGrid, -flip(yBaseGrid)])
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
    semilogy(optimValues.iteration,res,'ok')
    hold on
    xyLabelTex('\mathrm{iteration}','||R||_\infty')
    grid on
end

stop = 0;