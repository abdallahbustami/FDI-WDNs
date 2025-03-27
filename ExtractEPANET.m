%%ExtractEPANET.m
node_id = G.getNodeNameID;
PipeIndex = G.getLinkPipeIndex;
PumpIndex = G.getLinkPumpIndex;
JunctionIndex = G.getNodeJunctionIndex;
ReservoirIndex = G.getNodeReservoirIndex;
TankIndex = G.getNodeTankIndex;
LinkFromTo=G.getLinkNodesIndex;

% Run and extract hydraulics for verification (only fesibility)
H = G.getComputedHydraulicTimeSeries;
Flow = H.Flow;
Head = H.Head;
Pressure=H.Pressure;
Demand=H.Demand;
HeadLoss=H.HeadLoss;
JunctionCount = G.getNodeJunctionCount;
PipeCount = G.getLinkPipeCount;
PumpCount = G.getLinkPumpCount;
TankCount = G.getNodeTankCount;
ReservoirCount = G.getNodeReservoirCount;

NodesID=G.getNodeNameID;

%%Extracting data
%Hyadraulics, initial conditions and parameters from EPANET
PipeVelocity = H.Velocity(:, PipeIndex);
PipeLinkLength = G.getLinkLength(PipeIndex);
PipeLinkDiameter=G.getLinkDiameter(PipeIndex); %in inches
PipeLinkRoughnessCoeff=G.getLinkRoughnessCoeff(PipeIndex);
JunctionDemand = H.Demand(:, JunctionIndex); %demand at junction
JunctionElevation=G.getNodeElevations(JunctionIndex);
JunctionInitialHead=Head(1,JunctionIndex);
JunctionDemand=Demand(:,JunctionIndex);
PipeFlowRate = H.Flow(:, PipeIndex); %flow in pipe time series
PumpEnergy=H.Energy(PumpIndex);
if TankCount > 0
    TankElevation=G.getNodeElevations(TankIndex);
    TankInitialLevel=G.getNodeTankInitialLevel(TankIndex);
    TankMinLevel=G.getNodeTankMinimumWaterLevel(TankIndex);
    TankMaxLevel=G.getNodeTankMaximumWaterLevel(TankIndex);
    TankInitialHead=TankElevation+TankInitialLevel;
    TankMinHead=TankElevation+TankMinLevel;
    TankMaxHead=TankElevation+TankMaxLevel;
    TankDiameter=G.getNodeTankDiameter(TankIndex);
    TankArea=pi().*TankDiameter.*TankDiameter./4;
end
if ReservoirCount > 0
    ReservoirElevation=G.getNodeElevations(ReservoirIndex);
    ReservoirHead=Head(1,ReservoirIndex);
end

%% Pattern of price for pump 24 hrs
PatternIdx = find(strcmp(G.getPatternNameID, 'pr3'));
% PatternIdx = find(strcmp(d.getPatternNameID, 'pr4'));
PatternVal = G.getPattern;
PatternVal= PatternVal(PatternIdx,:);
PatternLength = G.getPatternLengths(PatternIdx);
PricePatternVal24h = repmat(PatternVal',24/PatternLength,1);

%% Get Demand 24 hrs
PatternIdx = find(strcmp(G.getPatternNameID, '1'));
% PatternIdx = find(strcmp(d.getPatternNameID, 'ModifiedPattern0'));
PatternVal = G.getPattern;
PatternVal= PatternVal(PatternIdx,:);
PatternLength = G.getPatternLengths(PatternIdx);
BaseDemand = G.getNodeBaseDemands{1};
BaseDemand = BaseDemand(JunctionIndex);
Demand72 = [];
DemandPatternVal24h = [];
JunctionDemand24=[];
for i = 1:PatternLength
    pat = repmat(PatternVal(i),24/PatternLength,1);
    DemandPatternVal24h = [DemandPatternVal24h; pat];
end

JunctionDemand24=[];
for JJ=1:JunctionCount
    JunctionDemand24temp=[];
    for i = 1 : length(JunctionDemand)
        if mod(i,60)==0
            JunctionDemand24temp=[JunctionDemand24temp; JunctionDemand(i,JJ)];
        end
    end
    JunctionDemand24=[JunctionDemand24, JunctionDemand24temp];
end

