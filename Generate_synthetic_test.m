clear all; clc; close all;

%%%%%%%%%%%%% settings %%%%%%%%%%%%%%%%%%
%rng(10)

% data settings
nrtrue =        0;
nrpol =         10;      % number of direct age constrints 
errorpol=       0.05;     % error (1 sigma) on measured age (%/100)
errorBe=        0.04;       % error (1 sigma) on measured 10Be (%/100)
errorAl=        0.04;      % error (1 sigma) on measured 10Be (%/100)
% OBS 3Kyr added error (1 sigma) on measured Be10 and Al26 (%/100)

n2= 2;                  % if n2=1 complex preburial history , else simple (norm distributed steady erotion rate)
pf=2;                   % prior distibution of erotion rate, if pf=1 then normal distribution, if pf=2 then uniform distribution. Only have effect if n2=2
simple=[0.001 0.03];    % if pf=1 then simple=[mean std], if pf=2 then simple=[min max] erotion rate (units?). Only have effect if n2=2

% define top and bottom depth and top age constraint
topd= 10;                  % Depth at the top of the age-depth model (m)
Age_at_top_depth = [10 50]; % (ka) Age at the top depth, if length(Age_at_top_depth)==1 then the age is locked at this age e.g. topd=0 and Age_at_top_depth=[0]
endd= 100;                 % Depth at the bottom of the age-depth model (m)

% Accumulation rate process 
meanacrate =    10;     % (kyr/m)
varac =         100;    % (kyr/m)
R =             [10];      % Correlation range (m) if length(ran)=1 fixed range, if length(ran)=2 then unifom distributed correlation range with ran=[min max]
restrue =       R(1)/10;     % Resolution, should be at least rantrue(1)/3 (m)

% Difine hiatus
Depth_of_hiatus = [35 55 75];          % (m) % the first defines the top/starting depth, if Depth_of_hiatus = [] no hiatus.
duration_of_hiatus = [10 500];         % Duration of the hiatus (ka)
h_correlated = 1;                      % if h_correlated=1 then the accumulations rates before and after each hiati are correlated, elles they are not. Note that the code does not work when h_correlated = 0, while including hiasuses and a variable R

% define constants
rho  =2;                               % density at the burial site (g/cm^3)
deBe =log(2)/1387000;                  % decay constant for 10Be (1/yr)
deAl =log(2)/705000;                   % decay constant for 26Al (1/yr)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% code starts from here

n =  1;                  % always keep as 1

% generate values
if nrtrue>0
    truedepth= sort(rand(nrtrue,1)*(endd-topd)+topd); % choose measured depth randomly
else 
    truedepth=[];
end
if nrpol>0
    poldepth= sort(rand(nrpol,1)*(endd-topd)+topd); % choose measured depth randomly
else 
    poldepth=[];
end

load('b.mat'); % load production rates



% cosmo constants
rho=2;                                % rock density at the burial site (g/cm^3)
deBe=log(2)/1387000;                  % decay constant for 10Be (1/yr)
deAl=log(2)/705000;                   % decay constant for 26Al (1/yr)

% data option
n = 1;                  % if n=1 fit both Al and Be, if n=2 Al/Be


% SETUP PRIOR
b2=varac/meanacrate; a2=meanacrate/b2;               % Gamma distiburtion scale and shape parameter
% depth vector
 dd=topd;lh=0;
if Depth_of_hiatus>0
    for i = 1:length(Depth_of_hiatus)
        dd=[dd dd(end):restrue:Depth_of_hiatus(i)];
        if dd(end)<Depth_of_hiatus(i)
            dd=[dd Depth_of_hiatus(i)];
        end
        lh(i)=length(dd);
    end

end
    dd=[dd dd(end):restrue:endd]; 
    if dd(end)<endd
        dd=[dd endd];
    end
    dd(1)=0;
lh(length(Depth_of_hiatus)+1)=length(dd); % index of hiatus
h_index=diff([0 lh-1])-1;

%sippi_set_path
 if h_correlated==1;
     forward.n_acr=length(dd)-length(Depth_of_hiatus)-2;
 else
forward.n_acr=h_index; % number of accumulation rates between each section with hiatus inbetween
 end
forward.n_cosmo = length(truedepth);
forward.truedepth=truedepth;
forward.poldepth=poldepth;
forward.ds=topd;
forward.production=b;
forward.data=n;
forward.preburial=n2;
forward.ran=R;
forward.rho=rho;
forward.debe=deBe;
forward.deal=deAl;
forward.dd=dd;
forward.topage=Age_at_top_depth;
forward.hc=h_correlated;
forward.h_index=h_index;

% SETUP PRIOR
b2=varac/meanacrate; a2=meanacrate/b2;               % Gamma distiburtion scale and shape parameter

clear prior
im=0;
if h_correlated==1;
    k2=1;
else
    k2=length(Depth_of_hiatus)+1;
end

for i=1:k2
im=im+1;
if length(R)<2 & R==0
    prior{im}.name=sprintf('Accummulation rates %d',im);
    prior{im}.type='gamma';
    prior{im}.x=1:forward.n_acr(i);
    prior{im}.a=ones(forward.n_acr(i),1).*a2;
    prior{im}.b=ones(forward.n_acr(i),1).*b2;
    %prior{im}.a(1)=50;
elseif length(R)<2 & R>0
    d_target = gaminv(rand(1,1000),a2,b2);
    prior{im}.name=sprintf('Accummulation rates %d',im);
    prior{im}.type='cholesky';
    prior{im}.x=[1:forward.n_acr(i)];
    prior{im}.Cm=sprintf('0.00001 Nug(0) + 0.99999 Gau(%g)',R/restrue);
    prior{im}.d_target = d_target;
elseif length(R)>1
    d_target = gaminv(rand(1,1000),a2,b2);
    prior{im}.name=sprintf('Accummulation rates %d',im);
    prior{im}.type='FFTMA';
    prior{im}.x=[1:forward.n_acr(i)]*restrue;
    prior{im}.Cm='1 Gau(10)';
    prior{im}.d_target = d_target;
end
end

if forward.n_cosmo>0
if forward.preburial==1 % Complex Pre-burial
    load BAl_fin;load NBe_fin;
    X=log(NBe_fin);xl='log(Be)';Y=log(NAl_fin./NBe_fin);yl='log(Al/Be)';
    NX=201;NY=201;[C,x_arr1,y_arr1]=histcounts2(X(:),Y(:),[NX,NY]);
    x_arr=(x_arr1(1:end-1)+x_arr1(2:end))/2;y_arr=(y_arr1(1:end-1)+y_arr1(2:end))/2;
    P=C./sum(C(:)); % this is the Complex 2D Pre-burial PDF

    for i=1:forward.n_cosmo;
        im=i+1;
        prior{im}.type='pdf2';
        prior{im}.pdf=P';
        prior{im}.pdf_x=x_arr;
        prior{im}.pdf_y=y_arr;
    end
    prior{im}.o_nscore.style='nearest';
else % Simple pre-burial history (steady erotion)
    im=im+1;
    if pf<2 % normal distribution
    prior{im}.x=1:forward.n_cosmo;
    prior{im}.name=sprintf('Erosion Rate');
    prior{im}.type='cholesky';
    prior{im}.m0=simple(1);
    prior{im}.Cm=sprintf('%f Nug(0)',simple(2).^2); 
    else % uniform distribution
    prior{im}.x=1:forward.n_cosmo;  
    prior{im}.type='uniform';
    prior{im}.name=sprintf('Erosion Rate');
    prior{im}.min=simple(1); %
    prior{im}.max=simple(2);
    end
end
end

if length(Age_at_top_depth)>1 % Age at the top
    im=im+1; 
    prior{im}.type='uniform';
    prior{im}.name='Age at top depth';
    prior{im}.min=Age_at_top_depth(1); % in meters
    prior{im}.max=Age_at_top_depth(2);
end

if length(Depth_of_hiatus)>0
    im=im+1; % hiatus
    prior{im}.x=1:length(Depth_of_hiatus); 
    prior{im}.type='uniform';
    prior{im}.name='Hiatus duration';
    prior{im}.min=duration_of_hiatus(1); % in meters
    prior{im}.max=duration_of_hiatus(2);
end

if  length(R)>1 % Variable correlation range
    im=im+1; 
    prior{im}.type='uniform';
    prior{im}.name='R';
    prior{im}.min=R(1); % in meters
    prior{im}.max=R(2);
    prior{im}.prior_master=k2;
end


% SETUP FORWARD
forward.forward_function='sippi_forward_CosmoChron';

% Test if the setup works 
m=sippi_prior(prior);               % generate a realization from the prior
%m{end}=[40 200 400]';
d=sippi_forward(m,forward,prior);   % calculate the forward respose to this realization


% setup for plotting
acr=[]; % set up prior accumulation rate
if h_correlated==1
        l1=1;
    for i=1:length(forward.h_index);
        l2=l1+forward.h_index(i)-1;
        acr=[acr Inf*ones(length(m{1}(:,1)),1) m{1}(:,l1:l2)];
        l1=l2+1;
    end
            im=1;
else
for i=1:length(forward.n_acr)
    acr=[acr Inf*ones(length(m{1}(:,1)),1) m{i}];
end
            im=i;
end

if forward.n_cosmo>0 % setup prior pre-burial
    if forward.preburial==1;
        preburial=[m{im+1:length(forward.truedepth)+im+1};];
            im=im+length(length(forward.n_acr));
    else
        preburial2=[m{im+1}];
        preburial=preburial2(:)';
            im=im+1;
    end
else
   preburial=[];
end

if length(forward.topage)==2 % setup top age
   top_age=m{im+1};
        im=im+1;
else
   top_age=forward.topage;
end

nrhiatus=sum(diff(forward.dd(2:end))==0);
if nrhiatus>0 % setup hiatus
    hiatus=m{im+1};
else
    hiatus=[];
end


tat = diff(dd).*acr; %  durating of each interval (kyr)
tat(1)=top_age;
if length(hiatus)>0
    tat(isnan(tat(1,:))) = hiatus(:);
end
trueage=[cumsum(tat')]; % true age (kyr)
depth=dd(2:end);


if n==1 & nrtrue>0
trueBe=d{1};   % true Be10 conten (Atons/g)
trueAl=d{2};   % true Al26 conten (Atons/g)
truestd = sqrt((1./trueBe).^2.*(trueAl*errorAl).^2 + (trueAl./trueBe.^2).^2.*(trueBe*errorBe).^2);    % true std on Al26/Be10 via error propagation
measuredBe=trueBe+randn(nrtrue,1).*trueBe*errorBe;   % true Be10 conten (Atons/g)
measuredAl=trueAl+randn(nrtrue,1).*trueAl*errorAl;   % true Al26 conten (Atons/g)
else
truedepth=[];trueBe=[];trueAl=[];truestd=[];measuredBe=[];measuredAl=[];eta=[];
end

if nrpol>0 & nrtrue>0
    errorpoli   = d{3}*errorpol;
    measuredpol= d{3}+randn(nrpol,1).*errorpoli; % measured age
elseif nrpol>0 & nrtrue==0
    errorpoli   = d{1}*errorpol;
    measuredpol= d{1}+randn(nrpol,1).*errorpoli; % measured age
else
     errorpoli   = [];
    measuredpol= [];
end

if nrtrue>0
if n2==2
 InitialBet=[b(1)./(deBe - rho*m{2}*b(5)) + b(2)./(deBe - rho*m{2}*b(6)) + b(3)./(deBe - rho*m{2}*b(7)) + b(4)./(deBe - rho*m{2}*b(8))];
 InitialAlt=[b(9)./(deAl - rho*m{2}*b(13)) + b(10)./(deAl - rho*m{2}*b(14)) + b(11)./(deAl - rho*m{2}*b(15)) + b(12)./(deAl - rho*m{2}*b(16))]; % initial 26Al and 10Be concentrations from constant erotion (Atoms/g), (eq. 3 in Ivy-OchsAndKober2008 with t=inf)
 eta=m{2};
else
   mvmvm=[m{2:1+length(forward.truedepth)};];
   InitialBet=exp(mvmvm(1,:));
   InitialAlt=exp(mvmvm(2,:)).*InitialBet;
   eta=[];
end
else
    InitialBet=[];
    InitialAlt=[];
    eta=[];
end

acrtrue=acr;
nrstepstrue=length(dd);

if forward.n_cosmo>0 & nrpol>0
   polage=d{3};
elseif forward.n_cosmo==0 & nrpol>0
   polage=d{1};
else
   polage=[];
end

if length(R)>1
truerage=m{im}; 
else
truerage=R;
end

truedd=dd;

save('TrueToFit.mat','truedd','trueage','truedepth','trueBe','trueAl','InitialBet','InitialAlt','truestd','measuredBe','measuredAl','errorBe','errorAl','acrtrue','endd','nrstepstrue','eta', 'polage','poldepth', 'measuredpol','errorpol','errorpoli','truerage')


% plot from here
if nrtrue>0
h(1)=subplot(121);hold on
errorbar(measuredAl./measuredBe,truedepth,[],[],truestd,truestd,'ok');
plot(trueAl./trueBe,truedepth,'sk','MarkerFaceColor','k');
xlabel('^2^6Al/^1^0Be');ylabel('Depth (m)'); ylim([0 endd]);set(gca,'Ydir','reverse');ylim([0 endd*1.05])
%plot(trueageattrue,truedepth,'sk','MarkerFaceColor','k');

h(2)=subplot(122); hold on
else
    figure;hold on;
end

% plot age-depth curve
nummmerer=1:length(dd);
hej=nummmerer(diff(dd(:))==0)-1;
shej=[1 hej+1];
ehej=[hej length(dd)-1];
for i = 1:length(shej)
  plot(trueage(shej(i):ehej(i)),dd(shej(i)+1:ehej(i)+1),'-k');hold on
end

% plot direct age constraints
if nrpol>0
errorbar(measuredpol,poldepth,[],[],errorpoli,errorpoli,'.k');
plot(polage,poldepth,'.k');end

% the simple burial ages using the the Metropolis algorithm
  intr=1e6; % number of iterations
 burnin=intr/5; % burnin period
for i=1:length(measuredAl);
[age2 eta2 Al2 Be2 accept likelihood]=SimpleBurialBalcoMCMC(intr,measuredAl(i),measuredAl(i)*errorAl,measuredBe(i),measuredBe(i)*errorBe, 2);
mcmcage=(median(age2(burnin:end))/1e3);
plus=prctile(age2(burnin:end)/1e3,84.1);
minus=prctile(age2(burnin:end)/1e3,15.9);
hold on
errorbar(mcmcage,truedepth(i),[],[],mcmcage-minus,plus-mcmcage,'.','Color',[1 1 1]*0.5)
end


xlabel('Age (kyr)');set(gca,'Xdir','reverse','Ydir','reverse');ylim([0 endd*1.05])
if nrtrue>0
set(h(1),'Position',[0.1    0.15    0.4    0.8],'fontsize',10);
set(h(2),'Position',[0.5    0.15    0.4    0.8],'fontsize',10,'yticklabel',[]);end
