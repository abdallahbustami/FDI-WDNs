%%PumpPowerConsumptionApprox.m
%%Simple optimization problem to obtain a convex approximation to the power
%%consumption function of pumps in the system
% clear
% close all

%%
%     ss=0.5; %pump speed
%pump charact
h0 = 393.7009;
alpha=3.7*10^(-6);
nu=2.59;
jj=1;

ii=50; %number of points
hPumpSeries =[];
qPumpSeries =[];
ssSeries=[];
ssmin=0.05;
ssmax=1;
ssstep=0.05;

for ss= ssmin : ssstep :  ssmax
    %get the max discharge q
    syms q0
    eqn = ss^2 * (h0 - alpha * (ss^(-1) * q0)^nu) == 0;
    q_Max(1,jj) = double(solve(eqn,q0,'Real',true));
    for i = 0 : ii
        ssSeries=[ssSeries , ss];
        qrand= i/ii*q_Max(1,jj);
        qPumpSeries=[qPumpSeries, qrand];
        hpump= ss^2 * (h0 - alpha * (ss^(-1) * qrand)^nu);
        hPumpSeries=[hPumpSeries, hpump];
    end
    jj=1+jj;
end

Nss=(ssmax-ssmin)/ssstep;
%%
%calculating the power for each of the points we have for discharge and
%head gain
SG=1; %Specific gravity for water
eta_Pump=0.85; %assume constant efficiency
for i = 1 : (ii+1)*(Nss+1)
    % PumpPowerkW(1,i) = qPumpSeries(1,i)*hPumpSeries(1,i)*SG/(3960*eta_Pump)*0.746; %in kW
    PumpPowerkW(1,i) = (qPumpSeries(1,i)/ConfigurationConstants.GPMperCFS)*hPumpSeries(1,i)*SG/(3960*eta_Pump)*0.746; %in kW

end

%%
cvx_begin
PowerSum=0;
clear PowerCoeff
clear PowerApp
variables PowerCoeff(6);
for i = 1 : (ii+1)*(Nss+1)
%     PowerApp(1,i) = ((PowerCoeff(1)*qPumpSeries(1,i)+PowerCoeff(2)*qPumpSeries(1,i)^2+PowerCoeff(3)*ssSeries(1,i)^2+PowerCoeff(4)*ssSeries(1,i)+...
%         +PowerCoeff(5))); %in kW
    PowerApp(1,i) = ((PowerCoeff(1)*qPumpSeries(1,i)+PowerCoeff(2)*qPumpSeries(1,i)^2+PowerCoeff(3)*ssSeries(1,i)^2+PowerCoeff(4)*ssSeries(1,i)+...
        +PowerCoeff(5)*ssSeries(1,i)*qPumpSeries(1,i)+PowerCoeff(6))); %in kW
end
P_Hussien=[2*PowerCoeff(3), PowerCoeff(6);
    PowerCoeff(6), 2*PowerCoeff(5)];
minimize norm(PowerApp-PumpPowerkW)
subject to
P_Hussien==semidefinite(2)
PowerApp>=zeros(1,length(PowerApp));
% PowerCoeff(1)>=0;
% PowerCoeff(3)>=0;
cvx_end

% save('PumpPowerConsumpAppCoeffs.mat','PowerCoeff')