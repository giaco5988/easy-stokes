%check metrics correlating physical soace and curvilinear coordinate

function [xGrid,yGrid] = checkMetricCurvilinearCheb(x,y,toleranceRemesh)

    % 
    mapping = @(t) t;
    t = linspace(0,1,numel(x))';

    % get chebfun
    x = chebfun(x,[0 1],'equi');
    y = chebfun(y,[0 1],'equi');

    %compute arclenght and metrics parameter
    l = computeTotalArcLengthSpectralCheb(x,y);
    xp = diff(x);    yp = diff(y);
    dl = sqrt(xp.^2+yp.^2)/l;

    % check if mesh position is different from physicsl one
    R = dl-diff(chebfun(mapping(t),[0 1],'equi'));
    doRemesh = sum(R.^2);

    % if the two are too different, flag to activate remesh
    if doRemesh>toleranceRemesh
        disp('Do remesh spectral for equi chebfun')
        options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','iter','MaxIter',1000);
        nonLinear = @(dof) remeshSolveCheb(dof,t,x,y,mapping);
        [tNew,~,~,output] = fsolve(nonLinear,t(2:end-1),options);
        if strcmp(output.message(1:17),'No solution found')
            error('remesh fails')
        end
        xGrid = x([0; tNew; 1]);
        yGrid = y([0; tNew; 1]);
    else
        xGrid = x(t);
        yGrid = y(t);
    end

end