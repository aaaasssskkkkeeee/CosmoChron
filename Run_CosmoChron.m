clear all;close all;clc;

% Make sure to download sippi https://github.com/cultpenguin/sippi and add to
% path by running: (in the command window):
% >> addpath sippi
% then run >> sippi_set_path

%% type-in data

% if you have no direct age constraints or comogenic nuclides data, then leave
% the constraints with square brackets

depth_of_direct_age_constraints = [];     % in meters (m)
direct_age_constraints = [];              % in thousand of years (ka)
uncertaintes_direct_age_constraints = []; % one sigma (ka)

depth_of_comogenic_nuclides = [];           % in meters
Al = [];                                    % Aluminium 26 (atom/g)
uncertaintes_26Al = [];                     % uncertainty on 26Al (one sigma) (atom/g)
Be = [];                                    % Berylium 10 (atom/g)
uncertaintes_10Be = [];                     % uncertainties on 10Be (one sigma) (atom/g)

%% settings

% file name
filename='CosmoChron_age-depth_model_1'

% pre-burial settings
n2= 1;                  % if n2=1 complex preburial history , if n2=2 simple (norm distributed steady erotion rate)
pf= 2;                  % prior distibution of erotion rate, if 1 normal, if 2 uniform, only have effect if n2=2
simple=[0.001 0.5];    % (cm/a) if pf=1 mean and std, if pf=2 min max erotion rate, only have effect if n2=2
% data option
n = 2;                  % if n=1 fit both Al and Be, if n=2 Al/Be

% define top and bottom depth and top age constraint
topd= 10;                  % Depth at the top of the age-depth model (m)
Age_at_top_depth = [10 20]; % age at the top (ka)
endd= 100;                 % Depth at the bottom of the age-depth model (m)

% accumulation process settings
meanacrate = 10;        % Mean inverse accumulation rate (kyr/m)
varac = 100;            % Variance on the inverse accumulation rate (kyr^2/m^2)
R =  [10];              % Correlation range (m) if length(R)=1 the fixed range, if length(R)=2 then unifom distributed correlation range with R=[min max] (variable R), R must be < endd-topd
res =  R(1)/6;          % resolution in depth, should be <=R(1)/3 (m)

% Difine hiatus
Depth_of_hiatus = [20 50 70];          % (m) % the first defines the top/starting depth, if Depth_of_hiatus = [] no hiatus.
duration_of_hiatus = [0 100];         % Duration of the hiatus (ka)
h_correlated = 1;                      % if h_correlated=1 then the accumulations rates before and after each hiati are correlated, elles they are not. obs if h_coorelation not 1 the R can not be a variable

% Extended Metropolish Sampling settings
options.mcmc.nite=1e4;               % Number if iterations
options.mcmc.i_sample=100;           % Number of interations between saving the sample
options.mcmc.i_plot=10000;           % Number of iterations between updating plots
options.mcmc.i_save_workspace=1e15;  % Number of iterations between saving the complete workspace
Anneal_length= 0.1;                  % Fraction of the sample length where the temperature annealation 
Adjust_steplength_length=0.11;       % Fraction of the sample length where the the steplength is adjustion,  Adjust_steplength_length must be grather than Anneal_length
% the burnin period is = Adjust_steplength_length*options.mcmc.nite by  default


% Cosmogenic nuclides constants
rho=2;                                % rock density at the burial site (g/cm^3)
deBe=log(2)/1387000;                  % decay constant for 10Be (1/yr)
deAl=log(2)/705000;                   % decay constant for 26Al (1/yr)
   

%% code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% code start from here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SETUP PRIOR
b2=varac/meanacrate; a2=meanacrate/b2;               % Gamma distiburtion scale and shape parameter
% depth vector
 dd=topd;lh=0;
if Depth_of_hiatus>0
    for i = 1:length(Depth_of_hiatus)
        dd=[dd dd(end):res:Depth_of_hiatus(i)];
        if dd(end)<Depth_of_hiatus(i)
            dd=[dd Depth_of_hiatus(i)];
        end
        lh(i)=length(dd);
    end
end
    dd=[dd dd(end):res:endd]; 
    if dd(end)<endd
        dd=[dd endd];
    end
    dd(1)=0;
lh(length(Depth_of_hiatus)+1)=length(dd);
h_index=diff([0 lh-1])-1;
load('b.mat'); % load production rates

%sippi_set_path
 if h_correlated==1;
     forward.n_acr=length(dd)-length(Depth_of_hiatus)-2;
 else
     forward.n_acr=h_index; % number of accumulation rates between each section with hiatus inbetween
 end
forward.n_cosmo = length(depth_of_comogenic_nuclides);
forward.truedepth=depth_of_comogenic_nuclides(:);
forward.poldepth=depth_of_direct_age_constraints(:);
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
    prior{im}.Cm=sprintf('0.00001 Nug(0) + 0.99999 Gau(%g)',R/res);
    prior{im}.d_target = d_target;
elseif length(R)>1
    d_target = gaminv(rand(1,1000),a2,b2);
    prior{im}.name=sprintf('Accummulation rates %d',im);
    prior{im}.type='FFTMA';
    prior{im}.x=[1:forward.n_acr(i)]*res;
    prior{im}.Cm='1 Gau(10)';
    prior{im}.d_target = d_target;
end
end

if forward.n_cosmo>0
if forward.preburial==1 % Complex Pre-burial
    load mc_source.mat;
    X=log(mc_source.NBe_fin);xl='log(Be)';Y=log(mc_source.NAl_fin./mc_source.NBe_fin);yl='log(Al/Be)';
    NX=201;NY=201;[C,x_arr1,y_arr1]=histcounts2(X(:),Y(:),[NX,NY]);
    x_arr=(x_arr1(1:end-1)+x_arr1(2:end))/2;y_arr=(y_arr1(1:end-1)+y_arr1(2:end))/2;
    P=C./sum(C(:)); % this is the Complex 2D Pre-burial PDF

    for i=1:forward.n_cosmo;
        im=im+1;
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
    prior{im}.prior_master=[1:length(Depth_of_hiatus)+1];
end

% SETUP DATA
nd=1; % index of data
if forward.n_cosmo>0
    if n == 1;
    data{nd}.d_std=uncertaintes_10Be(:);
    data{nd}.d_obs=Be(:);

    data{nd+1}.d_std=uncertaintes_26Al(:);
    data{nd+1}.d_obs=Al(:);
    else
    % Here I use error propagation
    data{nd}.d_std=sqrt((1./Be(:)).^2.*uncertaintes_26Al(:).^2 + (Al(:)./Be(:).^2).^2.*uncertaintes_10Be(:).^2);
    data{nd}.d_obs=Al(:)./Be(:);
    end
end

if length(forward.poldepth)>0
    if nd==1
        nd=0;end
    data{nd+1}.d_std=uncertaintes_direct_age_constraints(:);
    data{nd+1}.d_obs=direct_age_constraints(:);
end

% SETUP FORWARD
forward.forward_function='sippi_forward_CosmoChron';

%% Test if the setup works 
m=sippi_prior(prior);               % generate a realization from the prior
d=sippi_forward(m,forward,prior);   % calculate the forward respose to this realization
[logL]=sippi_likelihood(d,data)     % calculate the log likelihood

%%
% Run Extended Metropolis

% ANNEALING (TEMPERATURE AS A FUNCTION OF ITERATION NUMBER)
doAnneal=1;
if doAnneal==1;
    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=ceil(options.mcmc.nite*Anneal_length); %  iteration number when annealing stops
    options.mcmc.anneal.T_begin=50; % Start temperature for annealing
    options.mcmc.anneal.T_end=1; % End temperature for annealing
end

% Set step-length
for ip=1:length(prior)
    prior{ip}.seq_gibbs.step=1;
    prior{ip}.seq_gibbs.i_update_step_max=ceil(options.mcmc.nite*Adjust_steplength_length); %  iteration number when step size is locked
end

% Run Extended Metropolis sampling
options.txt=sprintf('D%d_PB%d',forward.data,forward.preburial);
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);

%%

%%%%%%%%%%%%%%%%%% plot the result %%%%%%%%%%%%%%%%%%%% 

nrburn=floor(prior{ip}.seq_gibbs.i_update_step_max/options.mcmc.i_sample)+1; % burn in period

for i=1:length(prior) % import data from file
[r{i},etype_mean{i},etype_var{i},reals_all{i},reals_ite{i}]=sippi_get_sample(options.txt,i);
end
 
% define parameters

nummmerera=1:length(dd); % for plotting hiatus
heja=nummmerera(diff(dd(:))==0)-1;
sheja=[1 heja+1];
eheja=[heja length(dd)-1];


acr=[]; % set up prior accumulation rate
if h_correlated==1
        l1=1;
    for i=1:length(forward.h_index);
        l2=l1+forward.h_index(i)-1;
        acr=[acr Inf*ones(length(reals_all{1}(:,1)),1) reals_all{1}(:,l1:l2)];
        l1=l2+1;
    end
            im=1;
else
for i=1:length(forward.n_acr)
    acr=[acr Inf*ones(length(reals_all{1}(:,1)),1) reals_all{i}];
end
            im=i;
end

if forward.n_cosmo>0 % setup prior pre-burial
    if forward.preburial==1;
        preburial=[reals_all{im+1:length(forward.truedepth)+im};];
            im=im+length(forward.truedepth);
    else
        preburial2=[reals_all{im+1}];
        preburial=preburial2(:)';
            im=im+1;
    end
else
   preburial=[];
end

if length(forward.topage)==2 % setup top age
   top_age=reals_all{im+1};
        im=im+1;
else
   top_age=forward.topage;
end

nrhiatus=sum(diff(forward.dd(2:end))==0);
if nrhiatus>0 % setup hiatus
    hiatus=reals_all{im+1};
else
    hiatus=[];
end


tat = diff(dd).*acr(nrburn:end,:); %  durating of each interval (kyr)
tat(:,1)=top_age(nrburn:end);
if length(hiatus)>0
    tat(:,isnan(tat(1,:))) = hiatus(nrburn:end,:);
    lmuh=prctile(hiatus(nrburn:end,:),15.9);muh=prctile(hiatus(nrburn:end,:),50);umuh=prctile(hiatus(nrburn:end,:),84.1);
    for i=1:length(Depth_of_hiatus)
            %subplot(1,length(Depth_of_hiatus),i) ; 
            figure; histogram(hiatus(nrburn:end,i),'Normalization','pdf','FaceColor',[1 1 1]*0.5);
            set(gca,'YTick',[]);xlabel('Duration (ka)');title(['Hiatus ' num2str(i) ': ' num2str(round(muh(i))) '+' num2str(round(umuh(i)-muh(i))) '-' num2str(round(muh(i)-lmuh(i))) ' ka'])
            xlim([duration_of_hiatus(1) duration_of_hiatus(2)])
    end
end
age=[cumsum(tat')]; % true age (kyr)
depth=dd(2:end);


% get forward response to data
clear albe polage
for i=1:length(reals_all{1}(:,1))
    if length(R)>1
        ost=1;else; ost=0; end
        for i2=1:length(reals_all)-ost
            mp(i2)={reals_all{i2}(i,:)};
        end
        if forward.preburial==1;
            for i3=1:length(forward.truedepth)
            mp(k2+i3)= {mp{k2+i3}'};end
        end

    d=sippi_forward(mp,forward);

    if forward.n_cosmo>0 & forward.data==1
        al(:,i)=d{2};
        be(:,i)=d{1};
        if length(forward.poldepth)>0
        polage(:,i)=d{3};   end
    elseif forward.n_cosmo==0 & length(forward.poldepth)>0 
        polage(:,i)=d{1};
    elseif forward.n_cosmo>0 & forward.data==2 & length(forward.poldepth)==0
        albe(:,i)=d{1};
    else
       albe(:,i)=d{1};
       polage(:,i)=d{2};
    end
end

if length(depth_of_comogenic_nuclides)>0
    if  n==1
        al=al(:,1:nrburn)';
        be=be(:,1:nrburn)';
    else
       albe=albe(:,1:nrburn)';
    end
end

% Age-depth plot
figure; hold on
if length(depth_of_comogenic_nuclides)>0
if n==1 % if complex pre-burial
h(1)=subplot(131);hold on
w=(endd-topd)/80;ggg=0.7;
for i=1:length(Al)
    fill([prctile(be(:,i),84.1),prctile(be(:,i),84.1), prctile(be(:,i),15.9), prctile(be(:,i),15.9)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w,depth_of_comogenic_nuclides(i)+w, depth_of_comogenic_nuclides(i)-w],[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.25)
    fill([prctile(be(:,i),97.7),prctile(be(:,i),97.7), prctile(be(:,i),2.3), prctile(be(:,i),2.3)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w,depth_of_comogenic_nuclides(i)+w, depth_of_comogenic_nuclides(i)-w],[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.15)
    plot([prctile(be(:,i),50) prctile(be(:,i),50)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'k','LineWidth',1.5)
    plot([prctile(be(:,i),84.1) prctile(be(:,i),84.1)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(be(:,i),15.9) prctile(be(:,i),15.9)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(be(:,i),97.7) prctile(be(:,i),97.7)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(be(:,i),2.3) prctile(be(:,i),2.3)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
        errorbar(Be,depth_of_comogenic_nuclides,[],[],uncertaintes_10Be,uncertaintes_10Be,'.k');
        xlabel('^1^0Be [atoms/g]');ylabel('Depth (m)'); ylim([0 max(depth)]);set(gca,'Ydir','reverse','XScale', 'log')
end
h(2)=subplot(132);hold on
w=(endd-topd)/80;ggg=0.7;
for i=1:length(Al)
    fill([prctile(al(:,i),84.1),prctile(al(:,i),84.1), prctile(al(:,i),15.9), prctile(al(:,i),15.9)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w,depth_of_comogenic_nuclides(i)+w, depth_of_comogenic_nuclides(i)-w],[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.25)
    fill([prctile(al(:,i),97.7),prctile(al(:,i),97.7), prctile(al(:,i),2.3), prctile(al(:,i),2.3)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w,depth_of_comogenic_nuclides(i)+w, depth_of_comogenic_nuclides(i)-w],[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.15)
    plot([prctile(al(:,i),50) prctile(al(:,i),50)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'k','LineWidth',1.5)
    plot([prctile(al(:,i),84.1) prctile(al(:,i),84.1)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(al(:,i),15.9) prctile(al(:,i),15.9)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(al(:,i),97.7) prctile(al(:,i),97.7)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(al(:,i),2.3) prctile(al(:,i),2.3)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
        errorbar(Al,depth_of_comogenic_nuclides,[],[],uncertaintes_26Al,uncertaintes_26Al,'.k');
        xlabel('^2^6Al [atoms/g]'); ylim([0 max(depth)]);set(gca,'Ydir','reverse','XScale', 'log')
end

else

h(1)=subplot(131);hold on
w=(endd-topd)/80;ggg=0.7;
for i=1:length(Al)
    fill([prctile(albe(:,i),84.1),prctile(albe(:,i),84.1), prctile(albe(:,i),15.9), prctile(albe(:,i),15.9)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w,depth_of_comogenic_nuclides(i)+w, depth_of_comogenic_nuclides(i)-w],[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.25)
    fill([prctile(albe(:,i),97.7),prctile(albe(:,i),97.7), prctile(albe(:,i),2.3), prctile(albe(:,i),2.3)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w,depth_of_comogenic_nuclides(i)+w, depth_of_comogenic_nuclides(i)-w],[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.15)
    plot([prctile(albe(:,i),50) prctile(albe(:,i),50)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'k','LineWidth',1.5)
    plot([prctile(albe(:,i),84.1) prctile(albe(:,i),84.1)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(albe(:,i),15.9) prctile(albe(:,i),15.9)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(albe(:,i),97.7) prctile(albe(:,i),97.7)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
    plot([prctile(albe(:,i),2.3) prctile(albe(:,i),2.3)],[depth_of_comogenic_nuclides(i)-w, depth_of_comogenic_nuclides(i)+w],'Color',[1 1 1]*ggg,'LineWidth',0.5)
end
ratioer=sqrt((1./Be(:)).^2.*uncertaintes_26Al.^2 + (Al(:)./Be(:).^2).^2.*uncertaintes_10Be(:).^2);
errorbar(Al./Be,depth_of_comogenic_nuclides,[],[],ratioer,ratioer,'.k');
xlabel('^2^6Al/^1^0Be');ylabel('Depth (m)'); ylim([0 max(depth)]);set(gca,'Ydir','reverse')
end

h(3)=subplot(133);hold on; 
for ia= 1: length(Depth_of_hiatus)+1
fill([prctile(age(sheja(ia):eheja(ia),:)',97.7), prctile(age(eheja(ia):-1:sheja(ia),:)',2.3)],[depth(sheja(ia):eheja(ia)), depth(eheja(ia):-1:sheja(ia))],[0 0.4470 0.7410],'EdgeColor',[1 1 1],'EdgeAlpha',0.1,'FaceAlpha',0.15)
fill([prctile(age(sheja(ia):eheja(ia),:)',84.1), prctile(age(eheja(ia):-1:sheja(ia),:)',15.9)],[depth(sheja(ia):eheja(ia)), depth(eheja(ia):-1:sheja(ia))],[0 0.4470 0.7410],'EdgeColor',[1 1 1],'EdgeAlpha',0.1,'FaceAlpha',0.25)
plot(prctile(age(sheja(ia):eheja(ia),:)',50),depth(sheja(ia):eheja(ia))','k','LineWidth',1.5);
end

% simple burial ages
intr=1e6; % number of iterations
burnin=intr/5; % define burnin period
for i=1:length(albe(1,:));
[age2 eta2 Al2 Be2 accept likelihood]=SimpleBurialBalcoMCMC(intr,Al(i),uncertaintes_26Al(i),Be(i),uncertaintes_10Be(i), 2);
mcmcage=(median(age2(burnin:end))/1e3);
plus=prctile(age2(burnin:end)/1e3,84.1);
minus=prctile(age2(burnin:end)/1e3,15.9);
hold on
errorbar(mcmcage,depth_of_comogenic_nuclides(i),[],[],mcmcage-minus,plus-mcmcage,'.','Color',[1 1 1]*0.5)
end

% direct age constraints
if length(depth_of_direct_age_constraints)>0
    errorbar(direct_age_constraints,depth_of_direct_age_constraints,[],[],uncertaintes_direct_age_constraints,uncertaintes_direct_age_constraints,'.k');
end

xlabel('Age (ka)');set(gca,'Xdir','reverse','Ydir','reverse')
ylim([0 max(depth)]);

if n==1
    set(h(1),'Position',[0.1    0.15    0.15    0.8],'fontsize',10);
    set(h(2),'Position',[0.25    0.15    0.15    0.8],'fontsize',10,'yticklabel',[]);
    set(h(3),'Position',[0.4    0.15    0.5    0.8],'fontsize',10,'yticklabel',[]); 
else
    set(h(1),'Position',[0.1    0.15    0.2    0.8],'fontsize',10);
    set(h(3),'Position',[0.3    0.15    0.6    0.8],'fontsize',10,'yticklabel',[]);
end


else
for i=1:length(eheja)
fill([prctile(age(sheja(i):eheja(i),:)',97.7), prctile(age(eheja(i):-1:sheja(i),:)',2.3)],[depth(sheja(i):eheja(i)), depth(eheja(i):-1:sheja(i))],[0 0.4470 0.7410],'EdgeColor',[1 1 1],'EdgeAlpha',0.1,'FaceAlpha',0.15)
fill([prctile(age(sheja(i):eheja(i),:)',84.1), prctile(age(eheja(i):-1:sheja(i),:)',15.9)],[depth(sheja(i):eheja(i)), depth(eheja(i):-1:sheja(i))],[0 0.4470 0.7410],'EdgeColor',[1 1 1],'EdgeAlpha',0.1,'FaceAlpha',0.25)
plot(prctile(age(sheja(i):eheja(i),:)',50),depth(sheja(i):eheja(i))','k','LineWidth',1.5);
end
if length(depth_of_direct_age_constraints)>0
    errorbar(direct_age_constraints,depth_of_direct_age_constraints,[],[],uncertaintes_direct_age_constraints,uncertaintes_direct_age_constraints,'.k');
end
ylim([0 endd*1.05])
xlabel('Age (ka)');ylabel('Depth (m)');set(gca,'Xdir','reverse','Ydir','reverse')
end
  

%% save the age-depth model

CosmoChron_95Lower=prctile(age',2.3)';
CosmoChron_68Lower=prctile(age',15.9)';
CosmoChron_Median=prctile(age',50)';
CosmoChron_68Upper=prctile(age',84.1)';
CosmoChron_95Upper=prctile(age',97.7)';
CosmoChron_Depth=depth';
T = table(CosmoChron_Depth,CosmoChron_95Lower,CosmoChron_68Lower,CosmoChron_Median,CosmoChron_68Upper,CosmoChron_95Upper);
writetable(T,filename)
