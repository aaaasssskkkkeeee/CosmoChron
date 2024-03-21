function [t_burial eta Al Be accept like_cur]=SimpleBurialBalcoMCMC(itr,measuredAl,errorAl,measuredBe,errorBe, rho)

%SimpleBurialBalcoMCMC(measuredAl(1),measuredAl(1)*errorAl(1),measuredBe(1),measuredBe(1)*errorBe(1), 2)
%A=griddedInterpolant(0:0.5:L,trueage)
%A(truedepth(1)) % in kyr

load('b.mat')
deAl=log(2)/705000; % (1/yr)
deBe=log(2)/1387000;
    Be=zeros(itr,1);
    Al=zeros(itr,1);
    like_cur=zeros(itr,1);
    t_burial=zeros(itr,1);
    eta=zeros(itr,1);

% initial
% t_burial(1)=Iage; % yr
% eta(1)=Ieta; %cm/yr
t_burial(1)=5e5; % yr
eta(1)=0.005; %cm/yr
like_cur(1)=-Inf;
ac=1;
RBe=[b(1)./(deBe - rho*eta(1)*b(5)) + b(2)./(deBe - rho*eta(1)*b(6)) + b(3)./(deBe - rho*eta(1)*b(7)) + b(4)./(deBe - rho*eta(1)*b(8))];
RAl=[b(9)./(deAl - rho*eta(1)*b(13)) + b(10)./(deAl - rho*eta(1)*b(14)) + b(11)./(deAl - rho*eta(1)*b(15)) + b(12)./(deAl - rho*eta(1)*b(16))];
Be(1)=RBe.*exp(-t_burial(1)*deBe);
Al(1)=RAl.*exp(-t_burial(1)*deAl);





for i=1:itr

%perturbate
t_burialt=t_burial(i)+randn*10e4;
etat=eta(i)+randn*0.001;
    
% forward
RBe=[b(1)./(deBe - rho*etat*b(5)) + b(2)./(deBe - rho*etat*b(6)) + b(3)./(deBe - rho*etat*b(7)) + b(4)./(deBe - rho*etat*b(8))];
RAl=[b(9)./(deAl - rho*etat*b(13)) + b(10)./(deAl - rho*etat*b(14)) + b(11)./(deAl - rho*etat*b(15)) + b(12)./(deAl - rho*etat*b(16))];
Bet=RBe.*exp(-t_burialt*deBe);
Alt=RAl.*exp(-t_burialt*deAl);

%likelihood
%if t_burialt>0
    liket=sum((-0.5*((Alt-(measuredAl)')./errorAl').^2)')+sum((-0.5*((Bet-(measuredBe)')./errorBe').^2)');
%else
%    liket=-Inf;
%end
%liket=normpdf(Alt,measuredAl,errorAl)*normpdf(Bet,measuredBe,errorBe);
%stdr=sqrt((1./measuredBe).^2.*(errorAl).^2 + (measuredAl./measuredBe.^2).^2.*(errorBe).^2); 
%liket=normpdf(Alt/Bet,measuredAl/measuredBe,stdr)*normpdf(Bet,measuredBe,errorBe);


%if liket/like_cur(i)>rand
if liket-like_cur(i)>log(rand)
    ac=ac+1;
    Be(i+1)=Bet;
    Al(i+1)=Alt;
    like_cur(i+1)=liket;
    t_burial(i+1)=t_burialt;
    eta(i+1)=etat;
else
    Be(i+1)=Be(i);
    Al(i+1)=Al(i);
    like_cur(i+1)=like_cur(i);
    t_burial(i+1)=t_burial(i);
    eta(i+1)=eta(i);   
end
accept=ac/itr;
end