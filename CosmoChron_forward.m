function [d1, d2, d3]=Cosmoforward(acr, preburial,top_age,hiatus,truedepth,poldepth,b,ds,n,n2, rho, deBe,deAl,dd,forward)


%%%%%%%%%%%%%%%%%%%%%%% simulated ages %%%%%%%%%%%%%%%%%%%%%%


ress=diff(dd); % depth reselution (m)
if ds>0
    acr(1)=top_age/dd(2);
end


if length(preburial)>0

nrtru=length(truedepth); % number of 26Al-10Be pairs 
acr2=0.1./acr;% accumulation rate (cm/yr)


initt1=(truedepth(end:-1:1)>dd)';
te2=sum(initt1(:,end:-1:1)); % interval including true age
missing=(truedepth' - dd(te2));

depthat=zeros(te2(end)+1,nrtru);%thic=zeros(te2(end),nrtru);%initt=zeros(te2(end)+2,nrtru); % Matrixes
for i=1:nrtru
     depthat((te2(end+1-i)):-1:1,i)= [0 cumsum(ress((te2(end+1-i)-1:-1:1)))]+missing(end+1-i);
end
thic=-diff(depthat);
initt=[ones(1,nrtru); depthat>0];

initt3=(initt-[initt(2:end,:); zeros(1,nrtru)]);

tat = thic(1:te2(end),:).*acr(1:te2(end))'*1000; %  durating of each interval (yr)
tat(1,:)=top_age*1000;

if length(hiatus)>0
idx = isnan(tat(:,1));
    if sum(idx)>0
    tat(idx,:) = hiatus(1:sum(idx)).*ones(sum(idx),nrtru)*1000; % add hiatus if any (yr)
    tat=tat.*initt(2:end-1,:);
    end
end


Pratesimat=zeros(te2(end)+1,nrtru,4);PrateAlsimat=zeros(te2(end)+1,nrtru,4);
for i=1:4
    Pratesimat(:,:,i)=      b(i)*exp(rho*100*depthat*b(i+4)).*initt(1:end-1,:);        % Balco, the 100 is to convert m to cm for 
    PrateAlsimat(:,:,i)=    b(i+8)*exp(rho*100*depthat*b(i+12)).*initt(1:end-1,:);
end

NBesimattru=zeros(1,nrtru);NAlsimattru=zeros(1,nrtru);
if n2==1 % complex preburial history ( preburial(2,:)=log(10Be), preburial(2,:)=log(26Al/10Be) ) 
    InitialBe = exp(preburial(1,:)); % initial 10Be (at/g)
    InitialAl = exp(preburial(2,:)).*InitialBe;% Initial 26Al
else % simple / stedy erotion pre burial history preburial = erootion rate (m/yr???)
    InitialBe=[b(1)./(deBe - rho*preburial*b(5)) + b(2)./(deBe - rho*preburial*b(6)) + b(3)./(deBe - rho*preburial*b(7)) + b(4)./(deBe - rho*preburial*b(8))];
    InitialAl=[b(9)./(deAl - rho*preburial*b(13)) + b(10)./(deAl - rho*preburial*b(14)) + b(11)./(deAl - rho*preburial*b(15)) + b(12)./(deAl - rho*preburial*b(16))]; % initial 26Al and 10Be concentrations from constant erotion (Atoms/g), (eq. 3 in Ivy-OchsAndKober2008 with t=inf)
end


for i=1:te2(end);
kk3=te2(end)+1-i;
NBesimattru=NBesimattru+initt3(kk3+1,:).*InitialBe; 
NBesimattru= (Pratesimat(kk3+1,:,1)./(deBe+rho*(acr2(kk3)')*b(5))) .* (exp(rho*(thic(kk3,:)*100)*b(5)) - exp(-deBe*tat(kk3,:)))... 
            +(Pratesimat(kk3+1,:,2)./(deBe+rho*(acr2(kk3)')*b(6))) .* (exp(rho*(thic(kk3,:)*100)*b(6)) - exp(-deBe*tat(kk3,:)))...
            +(Pratesimat(kk3+1,:,3)./(deBe+rho*(acr2(kk3)')*b(7))) .* (exp(rho*(thic(kk3,:)*100)*b(7)) - exp(-deBe*tat(kk3,:)))...
            +(Pratesimat(kk3+1,:,4)./(deBe+rho*(acr2(kk3)')*b(8))) .* (exp(rho*(thic(kk3,:)*100)*b(8)) - exp(-deBe*tat(kk3,:)))...
            + NBesimattru.*exp(-deBe*tat(kk3,:));

NAlsimattru=NAlsimattru+initt3(kk3+1,:).*InitialAl; 
NAlsimattru= (PrateAlsimat(kk3+1,:,1)./(deAl+rho*(acr2(kk3)')*b(13))) .* (exp(rho*(thic(kk3,:)*100)*b(13)) - exp(-deAl*tat(kk3,:)))...  
            +(PrateAlsimat(kk3+1,:,2)./(deAl+rho*(acr2(kk3)')*b(14))) .* (exp(rho*(thic(kk3,:)*100)*b(14)) - exp(-deAl*tat(kk3,:)))... 
            +(PrateAlsimat(kk3+1,:,3)./(deAl+rho*(acr2(kk3)')*b(15))) .* (exp(rho*(thic(kk3,:)*100)*b(15)) - exp(-deAl*tat(kk3,:)))... 
            +(PrateAlsimat(kk3+1,:,4)./(deAl+rho*(acr2(kk3)')*b(16))) .* (exp(rho*(thic(kk3,:)*100)*b(16)) - exp(-deAl*tat(kk3,:)))... 
            + NAlsimattru.*exp(-deAl*tat(kk3,:)); 

end


Beattrue=NBesimattru(:,end:-1:1)';
Alattrue=NAlsimattru(:,end:-1:1)';

if n== 2
    d1=Alattrue./Beattrue;
else
    d1=Beattrue;
    d2=Alattrue;
end
end

%%%%%%%%%%%%%%%%%%%%%% generate ages %%%%%%%%%%%%%%%%%%%%%%%

if length(poldepth)>0;

    duration = diff(dd(:)).*acr(:); %  durating of each interval (kyr)
    duration(1)=top_age;
    if length(hiatus)>0
        duration(isnan(duration)) = hiatus;
    end
    age=cumsum(duration)';  % sim age (kyr)
    %te1pol=(poldepth-hiatus_depth)/res;
    %te2pol=floor(te1pol);               % interval-1 including true age 
    %missingpol=(te1pol-te2pol)*res;     % depth span in te2+1 interval to true depth
    %polage=age(te2pol+1)+missingpol'.*acr(te2pol+2);
    %polage=interp1(dd,age,poldepth) 

    initt=(poldepth(end:-1:1)>dd)';
    te2pol=sum(initt(:,end:-1:1));
    missingpol=(poldepth'- dd(te2pol) );
    polage=age(te2pol-1) + missingpol.*acr(te2pol);

    if length(preburial)>0
        if n == 1
            d3=polage';
        else
            d2=polage';
        end
    else
        d1=polage';
    end
end

%end
