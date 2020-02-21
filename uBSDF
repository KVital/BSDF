clear all;close all; clc
%% Parameter section
MD=4;%years Mission duration (4 and 10 years)
Mrd=0.3;%m Mirror diameter
rho_target=2.53;%g/cm3 Zerodur density
CF=1;%Correction factor (shielding by structure). if CF==1 - no correction (worse scenario)
rho=2.5;%g/cm^3 constant density for all micrometeoroids
velo=20;%km/s impact velosity
alpha=0;%rad impact angle
lam=1.064e-6;%m wavelength 


%% Mictometeoroid amount
% Grun model
m=10.^(-(1:0.1:13));

F1=(2.2e+3*m.^0.306+15).^(-4.38);
F2=1.3e-9*(m+1e+11*m.^2+1e+27*m.^4).^(-0.36);
F3=1.3e-16*(m+1e+6*m.^2).^(-0.85);
Fs=3.1557e+7*(F1+F2+F3);%1/m2/year

Sm=pi*(Mrd/2)^2;% m^2 Mirror surface
Fs=Fs*MD*Mrd*Sm;
disp(['Number of impacts for  micrometeoriod mass >1e-12g  is ' num2str(round(Fs(find(m<=1e-12,1))))])
% figure(1)
% hold on
% plot(m,Fs,'-k','linew',2)
% grid on
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% set(gca,'fontsize',16)
% xlabel('Mass, [g]')
% ylabel('Number of M1 impacts')
% xlim([1e-13 1])

 

mass=m(2:end);%g fraction of the assumed number at a given mass should be assigned to the next lower mass (p3 SWALES)
NEI=(diff(Fs));% Number of expected impacts
NEI(NEI<=1/2.71828/2.71828)=0;% Cut on micrometeoroid size


Diam=(6*mass/(pi*rho)).^(1/3);%cm

figure(2)
hold on
plot(mass,NEI,'r','linew',2)
% plot(mass,diff(Fs),'-ok','linew',2)

grid on
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'fontsize',16)
xlabel('Mass, [g]')
ylabel('Number of impacts for a given mass')
xlim([ 1e-13 1e-6])
% legend('4 years','10 years')
% End of Stage 1
%% Crater Diameter

% James B. Heaney MLCD micrometeoroid optical damage analysis -- solar window 
D1=10^(1.569).*mass.^0.37;%cm
D2=10^(1.793).*mass.^0.396;%cm
D3=10^(1.485).*mass.^0.377;%cm


%% Particle/wall interaction model

% Gault Fechtig McHugh&Richardson Cour-Palais
K1=[1.08 6.0 1.28 1.06];
lambda=[1.071 1.13 1.2 1.06];
beta=[0.524 0.71 0 0.5];
gamma=[0.714 0.755 2/3 2/3];
xi=[0.714 0.755 2/3 2/3];
kappa=[-0.5 -0.5 0.5 0];

Kc=10;%Zerodur is brittle

for p=1:4
Dpw(p,:)=K1(p)*Kc*(Diam.^lambda(p))*(rho^beta(p))*(velo.^gamma(p))*(cos(alpha)^xi(p))*(rho_target^kappa(p));
end
% sum(Dpw.^2*pi.*repmat(NEI,[4 1]),2)
%% Uniter data for Stage 2
DCrater=[ D1; D2; D3; Dpw];

figure(3)
hold on
for cnt=1:7 
  plot(mass, DCrater(cnt,:),'linew',2)
end
grid on
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'fontsize',16)
xlabel('Mass, [g]')
ylabel('Crater diameter, [cm]')
legend('Horz 1','Horz 2','Horz 3','Gault','Fechtig','McHugh&R','Cour-Palais')
%% Ejected mass
% J-M Siguier J-C Mandeville. Test procedures to evaluate spacecraft
% materials ejecta upon hypervelocity impact

% Brittle target
K=ones(size(Diam));
K(Diam<0.001)=Diam(Diam<0.001)*1e+3;% comment this line for a ductile target

Ei=(mass/1000)*(velo*1000)^2/2;%Joule
Me=K.*7.41e-6.*sqrt(rho/rho_target).*Ei.^1.133.*cos(alpha)^2*1000;%g

betam=ones(size(Diam));
betam(Diam>0.01)=0.4;% comment this line for a ductile target
betam((Diam>0.0001)&(Diam<=0.01))=-0.3.*log10(Diam((Diam>0.0001)&(Diam<=0.01))/100)-0.8;% comment this line for a ductile target


Mcone=betam.*Me;
Mspalls=(1-betam).*Me;

figure(4)
hold on
plot(mass, Me./mass,'-k','linew',2)
plot(mass, Mcone./mass,'-r','linew',2)
plot(mass, Mspalls./mass,'-b','linew',2)

grid on

set(gca,'xscale','log')
 set(gca,'yscale','log')
set(gca,'fontsize',16)
xlabel('Mass \mumeteoroid, [g]')
ylabel('Ejected mass/Mass of \mumeteoroid')
legend('Total mass','Cone','Spalls','location','southeast')

%% Contamination model
Mt=sum(NEI.*Me);%g total emitted mass
VolEt=Mt*rho_target*(10000)^3; %um^3 volume of total emitted material

De=(6*mass/(pi*rho_target)).^(1/3)*10000;%um equivalent diameter

NEI01=(diff(Fs))/Sm*0.1; % particle distribution on 0.1m^2 surface
NEI01(NEI01<=1/2.71828/2.71828)=0;% Cut on micrometeoroid size

% max(Me(NEI01>0)) is the biggest ejected mass
CL=(6*max(Me(NEI01>0))/(pi*rho_target)).^(1/3)*10000;%um Cleanliness level is equal to the biggest particle diameter on the surface of 0.1m^2

% mass conservation
S=0.6:0.0001:1.5;
for cnt=1:length(S)

Nd=10.^(S(cnt).*(log10(CL).^2-log10(De).^2))/0.1*Sm;% MIL-STD-1246C particle distribution on the surface
fds=diff(Nd);

rez(cnt)=abs(sum(fds.*NEI(2:end).*((De(2:end)/2).^3))*4/3*pi-VolEt);
end
[~,psm]=min(rez);

disp(['Max ejected mass=' num2str(max(Me(NEI01>0))) 'g; S=' num2str(S(psm)) '; CL=' num2str(CL)])


%% BRDF Peterson model
thet=10.^(-(0:0.1:6))*pi/2;

BSDF=zeros([size(DCrater,1) length(thet)]);


figure(6)
hold on

for cnt=1:7

for cntx=1:length(NEI)%sum over all diameters
Nd=NEI(cntx)./Sm;%Density of the umeteoroids
d=DCrater(cnt,cntx)*1e-2;% D1 D2 D3 Dpw1 Dpw2 Dpw3 Dpw4
ld=(4/pi^4)^(1/3)*lam/d;
BSDF(cnt,:)=BSDF(cnt,:)+Nd*d^2/4*(1+(pi*d/(2*lam))^2*(1+(sin(thet)./ld).^2).^(-3/2));
end

plot(sin(thet),BSDF(cnt,:),'linew',2)
end

grid on
set(gca,'fontsize',16)
set(gca,'yscale','log')
 set(gca,'xscale','log')
xlabel('Sine of scatter angle')
ylabel('BSDF')
legend('Horz 1','Horz 2','Horz 3','Gault','Fechtig','McHugh&R','Cour-Palais')




%% ABg+R for optical simulation
 for p=1:7
bsdf=BSDF(p,:);
[xData, yData] = prepareCurveData( sin(thet), bsdf );
ft = fittype( [num2str(bsdf(end)) '*B/(B+x^2.55)+' num2str(bsdf(1))], 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'final';
opts.DiffMinChange=1e-14;
opts.TolX=1e-14;
opts.TolFun=1e-14;
opts.MaxIter=1000;
opts.StartPoint = [3e-9];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

B(p)=fitresult.B;
A(p)=bsdf(end)*B(p);
%  g=fitresult.g
Lmb(p)=bsdf(1);

 end
% 
  g=2.55;
R=Lmb*pi;% R is reflectance of the surface


format shortE
  disp('g=2.55') 
  disp('R   A   B ') 
for cnt=1:7
   disp([num2str(R(cnt))  '   '  num2str(A(cnt)) '   '  num2str(B(cnt)) ]) 
end



figure
hold on
for p=1:7
plot(thet,BSDF(p,:),'linew',2)
plot(thet,A(p)./(B(p)+sin(thet).^g)+Lmb(p),'--','linew',2)
end
grid on
set(gca,'fontsize',16)
set(gca,'yscale','log')
set(gca,'xscale','log')

%% BRDF Peterson model TIS
% thet1=0:0.000001:pi/2;
% 
% BSDF1=0;
% 
% 
% 
% for cnt=1:7
% 
% for cntx=1:length(NEI)%sum over all diameters
% Nd=NEI(cntx)./Sm;%Density of the umeteoroids
% d=DCrater(cnt,cntx)*1e-2;% D1 D2 D3 Dpw1 Dpw2 Dpw3 Dpw4
% ld=(4/pi^4)^(1/3)*lam/d;
% BSDF1=BSDF1+Nd*d^2/4*(1+(pi*d/(2*lam))^2*(1+(sin(thet1)./ld).^2).^(-3/2));
% end
% TIS(cnt)=sum(BSDF1.*cos(thet1).*sin(thet1))*2*pi*1e-6;
% end
% format shortE
% TIS
% 
% format shortE
% 
% for cnt=length(mass):-1:1
%    disp([num2str(mass(cnt))  '&'  num2str(Fs(cnt+1)) '&'  num2str(NEI(cnt)) '\\ \hline'  ]) 
% end
% 
%    disp([num2str(m(1))  '&'  num2str(Fs(1))  '&\\ \hline'  ]) 


