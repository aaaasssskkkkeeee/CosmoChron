clear all;close all;clc;


% define top and bottom depth and top age constraint
topd= 10;                  % Depth at the top of the age-depth model (m)
Age_at_top_depth = [10 60]; % (ka) Age at the top depth, if length(Age_at_top_depth)==1 then the age is locked at this age e.g. topd=0 and Age_at_top_depth=[0]
endd= 100;                 % Depth at the bottom of the age-depth model (m)

% accumulation process settings
meanacrate = 15;        % Mean inverse accumulation rate (kyr/m)
varac = 1500;            % Variance on the inverse accumulation rate (kyr^2/m^2)
R =  [5];              % Correlation range (m) if length(R)=1 the fixed range, if length(R)=2 then unifom distributed correlation range with R=[min max] (variable R), R must be < endd-topd
res =  R(1)/5;          % resolution in depth, should be <=R(1)/3 (m)

% Difine hiatus
Depth_of_hiatus = [];          % (m) Depth of the hiatuses. If Depth_of_hiatus = [] no hiatus.
duration_of_hiatus = [10 300];         % Duration of the hiatus (ka)
h_correlated = 1;   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% code start from here

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
forward.ds=topd;
forward.ran=R;
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


m=sippi_prior(prior);  % generate a relaization


% setup plotting
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
        acr=[acr Inf m{1}(l1:l2)];
        l1=l2+1;
    end
            im=1;
else
for i=1:length(forward.n_acr)
    acr=[acr Inf m{i}];
end
            im=i;
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
    tat(:,isnan(tat)) = hiatus;
end
age=[cumsum(tat')]; % true age (kyr)
depth=dd(2:end);


for i = 1:length(sheja)
  plot(age(sheja(i):eheja(i)),dd(sheja(i)+1:eheja(i)+1),'k','LineWidth',2);hold on
end
xlabel('Age (ka)');ylabel('Depth (m)');set(gca,'Xdir','reverse','Ydir','reverse')