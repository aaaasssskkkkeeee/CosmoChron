function [b]=AlBe_production_rate(elev,lat,lon)

% elev=5000;
% lat=47;
% lon=8;
clc
%nn = 100000;
%elev = 20000; %(m)
n_elev = 1;


%*********** Geographical settings ************
%lat = 47;                                                  % Latitude of sink site
%lon = 8;                                                % Longitude of sink site, negative if W
%sink_elev = 379;                                             % Elevation of sampling site [m asl.]
%max_elev = 4000;                                             % Maximum elevation of the possible source area [m asl.]

%************* Setup for calculations **************
%elev = linspace(sink_elev,max_elev,n_elev);                        % Elevation of possible source area [m asl.]
%eint_mc = logspace(log10(min_ero*1e-6),log10(max_ero*1e-6),nn);    % Erosion rates of source area [m/Myr]
%rho_pool = linspace(rho_min,rho_max,nn);                           % Rock density in source [g/cm3]
%rho_pool = 2; 
%pluck_depth = logspace(log10(pluck_min),log10(pluck_max),nn);      % Plucking depth pr. glaciation [m]

%elev_rough = linspace(min(elev),max(elev),n_elev);
elev_rough = elev; %(m)
Pspal_Be_rough = zeros(1,n_elev);
Pspal_Al_rough = zeros(1,n_elev);
%z = linspace(0,min([1.9e5 max_ero*1e-6*50*dt*100*max(rho_pool)]),1000);
z = linspace(0,20000,20000); % mass depth (g/cm^2) depth(m)*100(cm/m)*rho(g/cm^3)= mass depth for depth and rho

P_SLHL_Be = 4;                                               % 10Be production sea-level high-latitude [at/g/yr]
R_SLHL = 6.8;                                                % 26Al/10Be production ratio at SLHL
P_SLHL_Al = R_SLHL*P_SLHL_Be;                                % 26Al production SLHL
%dt = 2000;                                                   % Timesteps in calculations [yr]

%*********** Cosmo constants ************

Lspal_Be = 160;                                         % 10Be spallation attenuation length [g/cm^2]
Lspal_Al = 160;                                         % 26Al spallation attenuation length [g/cm^2]
Lnmc_Be = 1500;                                         % 10Be neutron-capture attenuation length [g/cm^2]
Lnmc_Al = 1500;                                         % 26Al neutron-capture attenuation length [g/cm^2]
Lfm_Be = 4320;                                          % 10Be fast muon attenuation length [g/cm^2]
Lfm_Al = 4320;                                          % 26Al fast muon attenuation length [g/cm^2]
%TBe = 1.387e6;                                          % 10Be half-life [yr]
%TAl = 0.705e6;                                          % 26Al half-life [yr]
%lambda_Be = log(2)/TBe;                                 % 10Be mean lifetime [yr^-1]
%lambda_Al = log(2)/TAl;                                 % 26Al mean lifetime [yr^-1]
gmr = -0.03417;                                         % Assorted constants for atmospheric pressure
dtdz = 0.0065;                                          % Lapse rate from standard atmosphere
c10.fC = 0.704;                                         % 10Be specific constants
c10.fD = 0.1828;
c10.Natoms = 2.006e22;
c10.fstar = 0.00157; % nyt fit fra Balco
c10.sigma190 = 37.8e-30; % nyt fit fra Balco
c26.fC = 0.296;                                         % 26Al specific constants
c26.fD = 0.6559;
c26.Natoms = 1.003e22;
c26.fstar = 0.0118;% nyt fit fra Balco
c26.sigma190 = 521e-30;% nyt fit fra Balco
                
%%%%%%%%%%%% Calcucation muon production rates %%%%%%%%%%%%%%

b_Be = zeros(6,n_elev);
b_Al = zeros(6,n_elev);
start_Be = [1 1 1 -1/Lnmc_Be -1/Lfm_Be -1/Lfm_Be]; % start til at fitte
start_Al = [1 1 1 -1/Lnmc_Al -1/Lfm_Al -1/Lfm_Al];
    

LSD_DG(lat, lon, elev_rough,0,1e7,-1,10);                % Calculating 10Be production scaling for spallation
    load LSDout;                                                % Obtaining 10Be production scaling for spallation
    P_time = P_SLHL_Be*mean(LSDout.Be);                         % Calculating 10Be production at height  for spallation
    Pspal_Be = P_time(1);                             
    
    LSD_DG(lat, lon, elev_rough,0,1e7,-1,26);                % Calculating 26Al productions scaling for spallation
    load LSDout;                                                % Obtaining 26Al production scaling for spallation
    P_time = P_SLHL_Al*mean(LSDout.Al);                         % Calculating 26Al production at height for spallation
    Pspal_Al = P_time(1);
    
    pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (elev*dtdz)) ) ); % Atmospheric pressure at elevation, Balco eller lifton
    
    P_mu_Be = P_mu_total_edited(z,pressure,c10,'no');           % Calculating 10Be muon production
    P_mu_Al = P_mu_total_edited(z,pressure,c26,'no');           % Calculating 26Al muon production
    
    objfcn = @(b,x) b(1)*exp(x*b(4)) + b(2)*exp(x*b(5)) + b(3)*exp(x*b(6)); % Multi-exponential function
    nlm_Be = fitnlm(z',P_mu_Be',objfcn,start_Be);               % Fitting multiexponential function to 10Be muon production
    b_Be(:,1) = nlm_Be.Coefficients.Estimate;
    start_Be = nlm_Be.Coefficients.Estimate; 
    nlm_Al = fitnlm(z',P_mu_Al',objfcn,start_Al);               % Fitting multiexponential function to 26Al muon production
    b_Al(:,1) = nlm_Al.Coefficients.Estimate;
    start_Al = nlm_Al.Coefficients.Estimate;


% subplot(121);plot(P_mu_Be,z,'O'); hold on
% plot(objfcn([ b_Be(:,1)],z),z,'LineWidth',2);
% set(gca,'Ydir','reverse'); xlabel('^1^0Be muon production rate (at/yr)'); ylabel(['Mass depth (g/cm^2)'])
% subplot(122);plot(P_mu_Al,z,'O'); hold on
% plot(objfcn([ b_Al(:,1)],z),z,'LineWidth',2);
%     set(gca,'Ydir','reverse'); xlabel(['^2^6Al m' ...
%         'uon production rate (at/yr)']);

    b=[Pspal_Be b_Be(1:3)' -1/Lspal_Be b_Be(4:6)' Pspal_Al  b_Al(1:3)' -1/Lspal_Be b_Al(4:6)'];
    save('b.mat','b')
end
