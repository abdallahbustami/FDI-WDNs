%%PipeHeadLoss_Linearization.m
function[QLinMidPoints,HLinMidPoints,PipeLinSegmSlope,PipeLinSegmIntcpt]=PipeHeadLoss_Linearization(PLength,PDiameter,PRoughnessCoeff,Npw,QLinMax,GPMperCFS)
%The Hazen-Williams equ. uses q in cfs for constant of 4.73 while EPANET
%gives flow in GPM
% GPM=CFS∗7.481∗60
HW_Const= 4.73;
QLinMax_CFS=QLinMax/GPMperCFS;
HLinMaxMin=HW_Const*PLength/((PRoughnessCoeff)^1.852)/(PDiameter^4.87)*QLinMax_CFS*abs(QLinMax_CFS^0.852);
QLinMidPoints=[0];
HLinMidPoints=[0];
PipeLinSegmSlope=[];
PipeLinSegmIntcpt=[];
for nn=1:(Npw/2)
    QLinMidPoint=nn/(Npw/2)*QLinMax;
    QLinMidPoint_CFS=QLinMidPoint/ConfigurationConstants.GPMperCFS;
    QLinMidPoints=[QLinMidPoints; QLinMidPoint];
    HLinMidPoint=HW_Const*PLength/((PRoughnessCoeff)^1.852)/(PDiameter^4.87)*QLinMidPoint_CFS*abs(QLinMidPoint_CFS^0.852);
    HLinMidPoints=[HLinMidPoints; HLinMidPoint];
    [PipeLinSegmCoeff]=polyfit([QLinMidPoint, QLinMidPoints(nn,1)],[HLinMidPoint,HLinMidPoints(nn,1)],1);
    PipeLinSegmSlope=[PipeLinSegmSlope; PipeLinSegmCoeff(1)];
    PipeLinSegmIntcpt=[PipeLinSegmIntcpt; PipeLinSegmCoeff(2)];
end
PipeLinSegmSlope=[flip(PipeLinSegmSlope); PipeLinSegmSlope];
PipeLinSegmIntcpt=[-flip(PipeLinSegmIntcpt); PipeLinSegmIntcpt];
QLinMidPoints=[-flip(QLinMidPoints(2:end)); QLinMidPoints];
HLinMidPoints=[-flip(HLinMidPoints(2:end)); HLinMidPoints];

% fprintf('Debug - Pipe Linearization:\n');
% fprintf('Length: %.2f ft\n', PLength);
% fprintf('Diameter: %.2f inches\n', PDiameter);
% fprintf('Roughness: %.2f\n', PRoughnessCoeff);
% fprintf('Max Flow (GPM): %.2f\n', QLinMax);
% fprintf('Max Head Loss: %.2f ft\n', HLinMaxMin);

end