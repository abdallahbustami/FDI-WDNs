%%PumpCurveApproximate.m
% H = a * Q^2 + b * Q + c * N^2 + d * N + e

%pump charact
h0 = 393.7009;
alpha=3.7*10^(-6);
nu=2.59;

ii=50; %number of points
HpumpAct =[];
qPumpSeries =[];
ssmin=0.05;
ssmax=1;
ssstep=0.05;
ssSeries=[];
jj=1;
for ss= ssmin : ssstep :  ssmax
    %get the max discharge q
    syms q0
    eqn = ss^2 * (h0 - alpha * (ss^(-1) * q0)^nu) == 0;
    q_Max(1,jj) = double(solve(eqn,q0,'Real',true));
    for i = 0 : ii
        ssSeries=[ssSeries, ss];
        qrand= i/ii*q_Max(1,jj);
        qPumpSeries=[qPumpSeries, qrand];
        hpump= ss^2 * (h0 - alpha * (ss^(-1) * qrand)^nu);
        HpumpAct=[HpumpAct, hpump];
    end
    jj=1+jj;
end

Nss=(ssmax-ssmin)/ssstep;
%
clear HpumpApp
clear PumpCoeff
cvx_begin
variables PumpCoeff(5);
for i = 1 : (ii+1)*(Nss+1)+1
    HpumpApp(1,i) = ((PumpCoeff(1)*qPumpSeries(1,i)^2+PumpCoeff(2)*qPumpSeries(1,i)+PumpCoeff(3)*ssSeries(1,i)^2+...
        +PumpCoeff(4)*ssSeries(1,i)+PumpCoeff(5))); %in ft
% HpumpApp(1,i) = ((PumpCoeff(1)*qPumpSeries(1,i)^2+PumpCoeff(2)*qPumpSeries(1,i)+PumpCoeff(3)*ssSeries(1,i)+...
%         +PumpCoeff(4))); %in ft
end
P_Hussien=[2*PumpCoeff(1), 0;
    0, 2*PumpCoeff(3)];
minimize norm(HpumpApp-HpumpAct)
subject to
HpumpApp>=zeros(1,length(HpumpApp));
P_Hussien==semidefinite(2);
PumpCoeff(1) >= 0;
cvx_end