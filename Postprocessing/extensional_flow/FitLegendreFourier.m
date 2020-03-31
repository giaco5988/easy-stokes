%fir dta to legendre polynomia

clear variables
close all

%physical parameters
Ca = 0.04;
lambda = 1;

%pathUpload = '~/Documents/C++/edge_state_extensional/resultsServer2/';
pathUpload = '~/Documents/C++/edge_state_extensional/resultsCaSmall/';
numberShapes = 1000;

%nodes
n = 200;

%options
Fourier = 1;                        % choose legendre or fourier
modes = 4;                          % choose number of modes
plotShape = 0;  PlotMesh = 0;       % plot shape option
PhaseSpace = 2;                     % choose phase space dimension

%colors
colors = 'brgkcm';
colors = repmat(colors,1,1000);

%parameters
IDshape = 1;
IDdelta = 1:16;
Round = 1;
V0 = 4/3*pi;
 
%loop on delta
for l = 1:numel(IDdelta)
    
    %initialize
    c1 = zeros(1,numberShapes);
    c2 = zeros(1,numberShapes);
    c3 = zeros(1,numberShapes);
    err = zeros(1,numberShapes);

%loop on shape
for k = 1:numberShapes
    
    display([num2str(k) ' of ' num2str(numberShapes)])
    
    %directory where to upload the results
    if numel(num2str(Ca))==4
        UploadResults = [pathUpload 'Ca=' num2str(Ca) '0000_lambda=' num2str(lambda) '.000000_IDshape=' num2str(IDshape) '_IDdelta=' num2str(IDdelta(l)) '_Round=' num2str(Round) '_elem=' num2str(n) '_dt=0.010000_loop=25000_RK=2_CPUs=1'];
    elseif numel(num2str(Ca))==3
        UploadResults = [pathUpload 'Ca=' num2str(Ca) '00000_lambda=' num2str(lambda) '.000000_IDshape=' num2str(IDshape) '_IDdelta=' num2str(IDdelta(l)) '_Round=' num2str(Round) '_elem=' num2str(n) '_dt=0.010000_loop=25000_RK=2_CPUs=1'];
    end
    name = [UploadResults '/drop' num2str(k-1) '.txt'];
    
    try
            A = importdata(name);
            YesBreak = 0;
    catch
           warning('No data')
           YesBreak = 1;
           break;
    end
    
    aIN = A.data(:,1)'; bIN = A.data(:,2)';
    
    %place in the center
    xcm = center_mass(aIN,bIN);
    aIN = aIN-xcm;
    a = aIN;    b = bIN;
    
    if plotShape==1
        
        figure(1)
        plot([a flip(a)],[b -flip(b)])
        if PlotMesh==1
            hold on
            plot([a flip(a)],[b -flip(b)],'o')
        end
        grid on
        axis equal
        axis([-3 3 -2 2])
        xlabel('x')
        ylabel('r')
        title(['IDshape=' num2str(IDshape) ' IDdelta=' num2str(IDdelta(l)) ' Round=' num2str(Round)])
        drawnow
        
    end
    
    r = sqrt(a.^2+b.^2);
    theta = atan(b./a);
    theta = theta + pi*(theta<0);

    %integration
    INT = ([diff(theta),0]+[0,diff(theta)])/2;

    %initialize
    f = zeros(1,modes);
    an= zeros(1,modes);
    bn= zeros(1,modes);

    if Fourier==2

        fitobject = fit(theta',r','poly3');

    else

        %compute modes
        SumP = 0;
        SumFourier = 0;
        for i = 1:modes

            if Fourier==0

                %polN = 2*i-1;
                polN = i-1;

                %legendre function
                %PPP = legendre(polN,cos(theta));
                PPP = legendre(polN,cos(theta));
                %legendre polinomials
                P = PPP(1,:);

                %compute coeff
                int = r.*P;
                f(i) = 2/pi*trapz(theta,int);

                if polN==0
                    SumP = f(i);
                else
                    SumP = SumP + f(i)*P;
                end

            elseif Fourier==1

                fourier = i-1;

                %take everything
                an(i) = 2/pi*trapz(theta,r.*cos(fourier*theta*2));
                bn(i) = 2/pi*trapz(theta,r.*sin(fourier*theta*2));

                if fourier==0
                    SumFourier = an(i)/2;
                else
                    SumFourier = SumFourier + an(i)*cos(fourier*theta*2) + bn(i)*sin(fourier*theta*2);
                end

            end

        end
        
        if Fourier==1;
            Vfourier = 2/3*pi*INT*(r'.^3.*sin(theta'));
            error = abs(Vfourier-V0)/V0;
        end
        err(k) = error;

    end
    
    %compute coefficients
%     cn = sqrt(an.^2+bn.^2);
%     %normalize
%     cn = cn/cn(1);
%     
%     col = [colors(l) 'x'];
%     
%     figure(2)
%     hold on
%     %plot(cn(2),cn(3),col)
%     plot3(cn(2),cn(3),cn(4),col)
%     hold off
%     grid on
%     xlabel('c_2')
%     ylabel('c_3')
%     drawnow
    
    %compute coefficients
    c1(k) = sqrt(an(2)^2+bn(2)^2)/sqrt(an(1)^2+bn(1)^2);
    c2(k) = sqrt(an(3)^2+bn(3)^2)/sqrt(an(1)^2+bn(1)^2);
    if modes==4
        c3(k) = sqrt(an(4)^2+bn(4)^2)/sqrt(an(1)^2+bn(1)^2);
    end

end
    
    if YesBreak==1
        c1 = c1(1:k-1);
        c2 = c2(1:k-1);
        if modes==4
            c3 = c3(1:k-1);
        end
        err = err(1:k-1);
    end    
        
    col = [colors(l) 'x-'];

    figure(2)
    if PhaseSpace==2
        plot(c1,c2,col)
        if k>1
            hold on
        end
        grid on
        title('Phase space')
        xlabel('c_1')
        ylabel('c_2')
        drawnow
    elseif PhaseSpace==3
        plot3(c1,c2,c3,col)
        if k>1
            hold on
        end
        grid on
        title('Phase space')
        xlabel('c_1')
        ylabel('c_2')
        zlabel('c_3')
        drawnow
    end
    
    %plot error
    figure(3)
    plot(err)
    if k>1
        hold on
    end
    grid on
    title('error truncation')
    xlabel('ite')
    ylabel('err_V')
    drawnow

end


%alpha = nthroot(4/3*pi/trapz(theta,2/3*pi*SumP.*sin(theta)),3);

%compare data to fit
% if Fourier==2
%     
%     %coeff
%     a = fitobject.p1;
%     b = fitobject.p2;
%     c = fitobject.p3;
%     d = fitobject.p4;
%     
%     f = a*theta.^3 + b*theta.^2 + c*theta.^1 + d;
%     
%     figure
%     plot(theta,r)
%     grid on
%     hold on
%     plot(theta,f)
%     xlabel('\theta')
%     ylabel('R')
%     legend('original','fitted','Location','Best')
%     
% else
% 
%     figure
%     plot(theta,r)
%     grid on
%     hold on
%     if Fourier==0
%         plot(theta,SumP)
%     elseif Fourier==1
%         plot(theta,SumFourier,'--')
%     end
%     xlabel('\theta')
%     ylabel('R')
%     legend('original','fitted','Location','Best')
% 
% end









