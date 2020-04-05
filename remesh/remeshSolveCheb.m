%non linear function for spectral remesh

function R = remeshSolveCheb(t,tOld,x,y,mapping)

    x = x([0; t; 1]);   y = y([0; t; 1]);

    % get chebfun
    x = chebfun(x,[0 1],'equi');
    y = chebfun(y,[0 1],'equi');

    %compute arclenght and metrics parameter
    l = computeTotalArcLengthSpectralCheb(x,y);
    xp = diff(x);    yp = diff(y);
    dl = sqrt(xp.^2+yp.^2)/l;

    % check if mesh position is different from physicsl one
    R = dl-diff(chebfun(mapping(tOld),[0 1],'equi'));
%     R = R(t);
    R = R(tOld);
    R(2:end-1);