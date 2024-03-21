% run
%AlBe_production_rate(1300,38.17,69.53); % elevation (m), latitude (^o), longitude (^o), values here: Kuldara
%AlBe_production_rate(4000,39.48,97.33); % tibet
AlBe_production_rate(379,47,8); % Mohlin

%%
clear all
load('b.mat')

rho=2;
z=linspace(0,100,100);

for i=1:4
    hold on
    Pratesimat(:,i) =   b(i)*exp(rho*100*z*b(i+4)); 
    PrateAlsimat(:,i)=  b(i+8)*exp(rho*100*z*b(i+12));
end

ro= 2.0;                              % dencity of rock (g/cm^2)
hat= [160 1500 4320];                 % attenuation length (spallation, negative nuons and fast muons) (g/cm^2)
prateAl26=[6.97*4.08-0.65 0.49 0.16]; % production rate of 26Al at z=0 from spallation, negative nuons and fast muons (atoms/(g*yr))
prateBe10=[4 0.06 0.02];              % production rate of 10Be (atoms/(g*yr))
deBe=log(2)/1387000;                  % decay constant for Be10 (1/yr)
deAl=log(2)/705000;                   % decay constant for Al26 (1/yr)
for i=1:3;
    PrateheisingerBe(:,i) = prateBe10(i).*exp(-ro*(z*100)/hat(i));
    PrateheisingerAl(:,i) = prateAl26(i).*exp(-ro*(z*100)/hat(i)); 
end

figure
plot(z,sum(PrateAlsimat')./sum(Pratesimat'));%hold on; plot(z,(Pratesimat(:,1)'));
hold on
plot(z,sum(PrateheisingerAl')./sum(PrateheisingerBe'));%hold on; plot(z,(PrateheisingerBe(:,1)'));
legend('Balco', 'Heisinger')
xlabel('Depth (m)'); ylabel('Production ^2^6Al/^1^0Be ratio')

figure
plot(z,sum(Pratesimat'));hold on; 
plot(z,(Pratesimat(:,1)'));
hold on
plot(z,sum(PrateheisingerBe'));
hold on; plot(z,(PrateheisingerBe(:,1)'));
legend('Balco', 'Balco without muon production','Heisinger','Heisinger without muon production')
xlabel('Depth (m)'); ylabel('Production rate (atoms/(g*y))')


figure
preburial=[0.0:0.00001:0.1]; rho=2;
LBe=[b(1)./(deBe - rho*preburial*b(5)) + b(2)./(deBe - rho*preburial*b(6)) + b(3)./(deBe - rho*preburial*b(7)) + b(4)./(deBe - rho*preburial*b(8))];
LAl=[b(9)./(deAl - rho*preburial*b(13)) + b(10)./(deAl - rho*preburial*b(14)) + b(11)./(deAl - rho*preburial*b(15)) + b(12)./(deAl - rho*preburial*b(16))];
loglog(LBe,(LAl./LBe),'k')
xlabel('^1^0Be (atoms/g)'); ylabel('^2^6Al/^1^0Be')
