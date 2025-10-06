function [d,forward,prior,data]=sippi_forward_CosmoChron(m,forward,prior,data,id,im)

if nargin<6,
    im=1;
end
if nargin<5,
    id=1;
end

acr=[]; % set up prior accumulation rate
if forward.hc==1
    l1=1;
    for i=1:length(forward.h_index)
        l2=l1+forward.h_index(i)-1;
        acr=[acr Inf m{1}(l1:l2)];
        l1=l2+1;
    end
else    
    for i=1:length(forward.n_acr)
        acr=[acr Inf m{i}];
    end    
end
im=length(forward.n_acr);

if forward.n_cosmo>0 % setup prior pre-burial
    if forward.preburial==1;
        preburial=[m{im+1:length(forward.truedepth)+im};];
            im=im+length(forward.truedepth);
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


if forward.n_cosmo>0 & length(forward.poldepth)>0 & forward.data==1
    [d{1},d{2},d{3}]=CosmoChron_forward(acr(:)', preburial,top_age,hiatus(:),forward.truedepth,forward.poldepth,forward.production,forward.ds,forward.data,forward.preburial,forward.rho,forward.debe,forward.deal,forward.dd);
elseif forward.n_cosmo==0 & length(forward.poldepth)>0 
    [d{1}]=CosmoChron_forward(acr(:)', preburial,top_age, hiatus(:),forward.truedepth,forward.poldepth,forward.production,forward.ds,forward.data,forward.preburial,forward.rho,forward.debe,forward.deal,forward.dd,forward);
elseif forward.n_cosmo>0 & forward.data==2 & length(forward.poldepth)==0
    [d{1}]=CosmoChron_forward(acr(:)', preburial,top_age,hiatus(:),forward.truedepth,forward.poldepth,forward.production,forward.ds,forward.data,forward.preburial,forward.rho,forward.debe,forward.deal,forward.dd);
else
    [d{1},d{2}]=CosmoChron_forward(acr(:)', preburial,top_age,hiatus(:),forward.truedepth,forward.poldepth,forward.production,forward.ds,forward.data,forward.preburial,forward.rho,forward.debe,forward.deal,forward.dd);
end
