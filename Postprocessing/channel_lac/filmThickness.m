%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Film thickness separating a bubble to the wall
%
% Gioele Balestra, LFMI / EPFL
%
% 6.10.2015
%
%   - start from a semicircle
%   - loop to fix L
%   - loop to fix Ca
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [x,y]=filmThickness(Caf,Lf)
clear variables; close all; clc;

%% Desired final parameters

Lf=1.9;
Caf=0.010;

%% Settings

FDorder=15; % order of finite difference scheme
NL=150;       % number of grid points for L-loop
deltaLStep=0.02;    % step for the change in L
deltaCaStep=0.005;   % step for the change in Ca
plotAll=false;    % plot everything
plotScalings=false;
plotVelocityProfile=true;
savedir='scalings';
mkdir(savedir);

%% Starting parameters, L-loop

L0=0.96;        % initial circle bubble radius
Ca0=0.08;         % capillary number

deltaL=sign(Lf-L0)*deltaLStep;     % continuation length for L-loop

% differentiation
[D1,D2,D3,D4,w,sL]=dif1D('fd',0,1,NL,FDorder);
Z=zeros(NL,NL); I=eye(NL);

%initial guess
sol=[L0/2*(1-cos(sL*pi)); 1/2-L0/2*sin(sL*pi); sL*pi-pi/2; L0/2*pi; L0];
Ca = Ca0;
dir=[zeros(3*NL+1,1);1];   % initial direction

disp('Continuation loop over L')
ind=0;quitcon=0;Qvar=0;count=0;
while ~quitcon%&ind<200
    
    if(Qvar>1e-6 || count>500) % interpolate to finer grid if flow rate not conserved or many iterations
        [sol,dir,D1,D2,D3,w,sL,Z,I,NL,FDorder] = interpolsol(sol,dir,sL,2*NL,FDorder);
    end
    
    solprevL=sol;
    sol=sol+dir*deltaL; % new prediction of solution
    
    % Newton iterations
    disp('Newton loop')
    quit=0;count=0;
    while ~quit
        
        % the present solution and its derivatives
        x=sol(1:NL);
        y=sol(NL+1:2*NL);
        th=sol(2*NL+1:3*NL);
        l=sol(3*NL+1);
        L=sol(3*NL+2);
        
        D1y = D1*y;
        
        D1th = D1*th;
        D2th = D2*th;
        D3th = D3*th;
        
        % nonlinear function
        f=[ ...
            l*cos(th)-D1*x; ...
            l*sin(th)-D1*y; ...
            1/(3*Ca)*y.^3/l^3.*(tan(th).*D1th.*D2th+D3th)+...
            1/Ca*y.^2/l^3.*D1y.*D2th-cos(th)/l.*D1y; ...
            x(NL)-L;...
            dir'*(sol-solprevL)-deltaL];
        
        % check if desired length obtained
        if(sign(Lf-L0)*(L-Lf)>=0)
            quitcon=1;
            LloopFinished=1;
            break;
        end
        
        
        % analytical jacobian
        A3y = diag(1/Ca*y.^2/l^3.*(tan(th).*D1th.*D2th+D3th))*I +...
            diag(1/Ca*2*y/l^3.*D1y.*D2th)*I +...
            diag(1/Ca*y.^2/l^3.*D2th)*D1+...
            diag(-cos(th)./l)*D1;
        
        A3th = diag(1/(3*Ca)*y.^3/l^3)*...
            (diag((1+tan(th).^2).*D1th.*D2th)*I +...
            diag(tan(th).*D2th)*D1 +...
            diag(tan(th).*D1th)*D2 + D3)+...
            diag(1/Ca*y.^2/l^3.*D1y)*D2+...
            diag(sin(th)./l.*D1y)*I;
        
        A3l = (-1/Ca*y.^3/l^4.*(tan(th).*D1th.*D2th+D3th))+...
            (-1/Ca*3*y.^2/l^4.*D1y.*D2th)+...
            (cos(th)./l^2.*D1y);
        
        
        A=[-D1,   Z,  -l*diag(sin(th)),   cos(th), Z(:,1); ...
            Z,   -D1, l*diag(cos(th)), sin(th), Z(:,1); ...
            Z,   A3y,  A3th, A3l, Z(:,1); ...
            I(NL,:),  Z(1,:), Z(1,:), 0, -1;
            dir'];
        
        
        % Boundary conditions
        loc = [1 NL+1 2*NL+1 3*NL 2*NL];
        f(loc)=[x(1)-0; th(1)+pi/2; y(1)-.5; y(NL)-.5; th(NL)-pi/2];
        A(loc,:)=[ ...
            I(1,:),Z(1,:),Z(1,:),0,0; ...
            Z(NL,:),Z(1,:),I(1,:),0,0; ...
            Z(NL,:),I(1,:),Z(1,:),0,0; ...
            Z(NL,:),I(NL,:),Z(1,:),0,0; ...
            Z(NL,:),Z(NL,:),I(NL,:),0,0];
        
        %        figure(4), imagesc(A), title('Jacobian')
        %        figure(5), plot(f), title('f')
        
        % convergence test
        %figure(6)
        %plot(x,y,'b-'); grid on; axis equal; drawnow;
        res=norm(f);
        disp([num2str(count) '  ' num2str(res)]);
        if count>1000||res>1e5; disp('no convergence');quitcon=1; break; end
        if res<1e-5; quit=1; disp('converged'); continue; end
        
        % Newton step
        sol=sol-A\f;
        count=count+1;
        %pause
    end
    
    % ploting the present solution
    if ~quitcon;
        Q = 1/Ca*1./(l^2.*cos(th)).*D2th.*(y.^3/3)-y; % flow rate at every x
        Qvar = std(Q(10:end-10))/abs(mean(Q(10:end-10))); % normalized std deviation of the flow rate
        
        figure(1)
        subplot(2,3,1:3); plot(x,y,'-k');
        xlabel('x'); ylabel('y');title('Shape of the interface')
        axis equal;
        subplot(2,3,4); plot(ind,L,'b.');
        xlabel('ind'); ylabel('L');
        hold on
        subplot(2,3,6); plot(ind,Qvar,'b.');
        xlabel('ind'); ylabel('std(Q)/|mean(Q)|');
        hold on;
        
        if plotAll
            figure(2)
            subplot(1,3,1); plot(1:1:NL,x);
            ylabel('x')
            subplot(1,3,2); plot(1:1:NL,y);
            ylabel('y')
            subplot(1,3,3); plot(1:1:NL,th);
            ylabel('th');
            hold off;

            figure(3)
            subplot(1,3,1);plot(1:1:NL,D1th);
            ylabel('D1*th')
            subplot(1,3,2);plot(1:1:NL,D2th);
            ylabel('D2*th')
            subplot(1,3,3);plot(1:1:NL,D3th);
            ylabel('D3*th');
            hold off;
        end
        drawnow;
    end
    
    % New direction
    dir=A\[zeros(3*NL+1,1);1]; % new direction
    dir=dir/norm(dir); % normalization
    
    ind=ind+1;
end
nfig = size(get(0,'Children'),1);

%% Starting parameters, Ca-loop

disp('Continuation loop over Ca')

NCa=NL;         % number of grid points
deltaCa=sign(Caf-Ca0)*deltaCaStep;    % continuation length for Ca-loop
deltaCaChanged=0;

dir=[zeros(3*NL+1,1);1];   % initial direction

[sol,dir,D1,D2,D3,D4,w,sCa,Z,I,NCa,FDorder] = interpolsol(solprevL,dir,sL,NCa,FDorder); % interpolate to finer grid if necessary (NCa~=NL)

sol(end)=Ca0; % change the last variable of sol to Ca0

ind=0;quitcon=0;Qvar=0;count=0;
while (~quitcon)&&LloopFinished%&ind<10000
    
    if(Qvar>1e-6 || count>500)
        [sol,dir,D1,D2,D3,w,sCa,Z,I,NCa,FDorder] = interpolsol(sol,dir,sCa,NCa+50,FDorder); % interpolate to finer grid if flow rate not conserved or many iterations
    end
    
    if (Ca<5*10^(-4) && ~deltaCaChanged)
        deltaCa = deltaCa/2; % reduce step size for small Ca
        deltaCaChanged=1;
        disp('NOTE: deltaCaChanged');
    end
    
    solprevCa=sol;
    sol=sol+dir*deltaCa; % new prediction of solution
    
    % Newton iterations
    disp('Newton loop')
    quit=0;count=0;
    while ~quit
        
        % the present solution and its derivatives
        x=sol(1:NCa);
        y=sol(NCa+1:2*NCa);
        th=sol(2*NCa+1:3*NCa);
        l=sol(3*NCa+1);
        Ca=sol(3*NCa+2);
        
        
        D1y = D1*y;
        
        D1th = D1*th;
        D2th = D2*th;
        D3th = D3*th;
        
        % nonlinear function
        f=[ ...
            l*cos(th)-D1*x; ...
            l*sin(th)-D1*y; ...
            1/(3*Ca)*y.^3/l^3.*(tan(th).*D1th.*D2th+D3th)+...
            1/Ca*y.^2/l^3.*D1y.*D2th-cos(th)/l.*D1y; ...
            x(NCa)-L;...
            dir'*(sol-solprevCa)-deltaCa];
        
        % check if desired Ca obtained
        if(sign(Caf-Ca0)*(Ca-Caf)>=0)
            quitcon=1;
            break;
        end
        
        % analytical jacobian
        A3y = diag(1/Ca*y.^2/l^3.*(tan(th).*D1th.*D2th+D3th))*I +...
            diag(1/Ca*2*y/l^3.*D1y.*D2th)*I +...
            diag(1/Ca*y.^2/l^3.*D2th)*D1+...
            diag(-cos(th)./l)*D1;
        
        A3th = diag(1/(3*Ca)*y.^3/l^3)*...
            (diag((1+tan(th).^2).*D1th.*D2th)*I +...
            diag(tan(th).*D2th)*D1 +...
            diag(tan(th).*D1th)*D2 + D3)+...
            diag(1/Ca*y.^2/l^3.*D1y)*D2+...
            diag(sin(th)./l.*D1y)*I;
        
        A3l = (-1/Ca*y.^3/l^4.*(tan(th).*D1th.*D2th+D3th))+...
            (-1/Ca*3*y.^2/l^4.*D1y.*D2th)+...
            (cos(th)./l^2.*D1y);
        
        A3Ca = (-1/(3*Ca^2).*y.^3./l^3.*(tan(th).*D1th.*D2th+D3th))+...
            (-1/Ca^2.*y.^2./l.^3.*D1y.*D2th);
        
        
        A=[-D1,   Z,  -l*diag(sin(th)),   cos(th), Z(:,1); ...
            Z,   -D1, l*diag(cos(th)), sin(th), Z(:,1); ...
            Z,   A3y,  A3th, A3l, A3Ca; ...
            I(NCa,:),  Z(1,:), Z(1,:), 0, 0;
            dir'];
        
        
        % Boundary conditions
        loc = [1 NCa+1 2*NCa+1 3*NCa 2*NCa];
        f(loc)=[x(1)-0; th(1)+pi/2; y(1)-.5; y(NCa)-.5; th(NCa)-pi/2];
        A(loc,:)=[ ...
            I(1,:),Z(1,:),Z(1,:),0,0; ...
            Z(NCa,:),Z(1,:),I(1,:),0,0; ...
            Z(NCa,:),I(1,:),Z(1,:),0,0; ...
            Z(NCa,:),I(NCa,:),Z(1,:),0,0; ...
            Z(NCa,:),Z(NCa,:),I(NCa,:),0,0];

%             figure(4), imagesc(A), title('Jacobian')
%             figure(5), plot(f), title('f')

        % convergence test
        %figure(6)
        %plot(x,y,'b-'); grid on; axis equal; drawnow;
        res=norm(f);
        disp([num2str(count) '  ' num2str(res)]);
        if count>1000||res>1e5; disp('no convergence');quitcon=1; break; end
        if res<1e-5; quit=1; disp('converged'); continue; end
        
        % Newton step
        sol=sol-A\f;
        count=count+1;
        %pause
    end
    
    % ploting the present solution
    if ~quitcon;
        Q = 1/Ca*1./(l^2.*cos(th)).*D2th.*(y.^3/3)-y;    % flow rate at every x location
        Qvar = std(Q(10:end-10))/abs(mean(Q(10:end-10)));    % normalized std deviation of the flow rate
        
        idxHinfAll = find((abs(D1th)<1e-4) & (abs(D2th)<1e-1));     % indeces of film with constant height
        [~,idxMin] = min(abs(D1th(idxHinfAll))+abs(D2th(idxHinfAll)));
        idxHinf = idxHinfAll(idxMin);
        
        Hinf = y(idxHinf);
        
        idxUc = 2;
        Uc = -1/Ca*1/(l^2*cos(th(idxUc)))*D2th(idxUc)*(y(idxUc).^2./2-y(idxUc).*y(idxUc))-1 +1; % centerline velocity in lab frame
        Umean = 2/3*(Uc); % mean velocity in lab frame
        
        Rf = 1/(2*(1+1.79*(3*Ca)^(2/3))); % asymptotic front radius
        Cx = x(end)-Rf;
        Cy = y(end);
        xf = Cx+Rf*cos(linspace(0,pi,100));
        yf = Cy-Rf*sin(linspace(0,pi,100));
        
        figure(1)
        hold on;
        subplot(2,3,1:3);plot(x,y,'-k',x(idxHinf),y(idxHinf),'g*');
        xlabel('x'); ylabel('y');title('Shape of the interface')
        axis equal;
        subplot(2,3,5); plot(ind,Ca,'g.');
        xlabel('ind'); ylabel('Ca');
        hold on;
        subplot(2,3,6); plot(ind,Qvar,'g.');
        xlabel('ind'); ylabel('std(Q)/|mean(Q)|');
        hold on;
        
        if plotScalings
            figure(nfig+1)
            loglog(Ca,Hinf,'ok','Markersize',5);
            xlabel('Ca'); ylabel('h_{inf}');
            title('Uniform film height')
            hold on; grid on;

            figure(nfig+2)
            loglog(Ca,(1-Umean)/1,'ok','Markersize',5);
            xlabel('Ca'); ylabel('(U-Umean)/U');
            title('Bubble velocity overshoot')
            hold on; grid on;

            figure(nfig+3)
            plot(x,y,'-k',xf,yf,'--r')
            xlabel('x'); ylabel('y');
            xlim([3/4*L L]);
            title('Front cap')
            axis equal;
            legend('numerics','asymptotics');
        end
        
        if plotAll
            figure(2)
            subplot(1,3,1); plot(1:1:NCa,x);
            ylabel('x')
            subplot(1,3,2); plot(1:1:NCa,y);
            ylabel('y')
            subplot(1,3,3); plot(1:1:NCa,th);
            ylabel('th');
            hold off;

            figure(3)
            subplot(1,3,1);plot(1:1:NCa,D1th);
            ylabel('D1*th')
            subplot(1,3,2);plot(1:1:NCa,D2th);
            ylabel('D2*th')
            subplot(1,3,3);plot(1:1:NCa,D3th);
            ylabel('D3*th');
            hold off;
        end
        drawnow;
    end
    
    % New direction
    dir=A\[zeros(3*NCa+1,1);1]; % new direction
    dir=dir/norm(dir); % normalization
    
    ind=ind+1;
end

%% Savings

save([savedir '/data.mat']);

figure(1)
print('-depsc','-loose','-r300',[savedir '/loopsHistory.eps'])

CaValues=linspace(Ca,Ca0,100);

if plotScalings
    figure(nfig+1)
    loglog(CaValues,0.5*1.3375.*CaValues.^(2/3),'r--',CaValues,0.5*1.3375.*CaValues.^(2/3)./(1+1.3375*2.5.*CaValues.^(2/3)),'g-','Linewidth',2);
    legend('numerics','1.3375*Ca^{2/3}','1.3375*Ca^{2/3}/(1+1.3375*2.5*Ca^{2/3})');
    print('-depsc','-loose','-r300',[savedir '/hInf.eps'])

    figure(nfig+2)
    loglog(CaValues,2*0.5*1.3375.*CaValues.^(2/3),'r--',CaValues,2*0.5*1.3375.*CaValues.^(2/3)./(1+1.3375*2.5.*CaValues.^(2/3)),'g-','Linewidth',2);
    legend('numerics','2*1.3375*Ca^{2/3}','2*1.3375*Ca^{2/3}/(1+1.3375*2.5*Ca^{2/3})');
    print('-depsc','-loose','-r300',[savedir '/Uovershoot.eps'])

    figure(nfig+3)
    print('-depsc','-loose','-r300',[savedir '/frontCap.eps'])
end

if plotVelocityProfile
    velocityProfile(x,y,th,l,Ca,L,D2th);
end

figure
plot(x,y)
axis equal

