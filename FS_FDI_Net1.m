%% Main code file for FS-FDI Implementation on Net1 (Single Pump Code)
clc; clear; close all;

%% Load EPANET Input file
G = epanet('Net1.inp');

%% Loading main important codes
% Extract results from EPANET and save EPANET file
ExtractEPANET

%constants between diff units
SecinMin = 60;
GPM2CFPS = 0.0022280092592593;

%Getting the system dimensions and elements index
%Components counts
nxH = JunctionCount + TankCount + ReservoirCount + PipeCount + PumpCount; %total number of elements
nLH=PipeCount + PumpCount;
nNH=JunctionCount + TankCount + ReservoirCount;

%Fitting Index for the problem
%We start by Nodes then Links, This makes it easier
IndexReservoir_H = ReservoirIndex;
IndexJunction_H = JunctionIndex;
IndexTank_H = TankIndex;
IndexPipe_H = PipeIndex+double(nNH);
if PumpCount > 0
    IndexPump_H = nxH;
end

%% Pipe linearization, pump curve approximation and pump power consumption
% Steps to perform as presteps of the optimization problem
% First we start with the pipe headloss curve piecewise linearization
% We need to specify how many pieces we want and for each linearization
% segment, the slope and the intercept along with the lines boundaries
Npw=4; %Choose an even number ALWAYS because the curve passes by (0,0)
QLinMax=1500; %in GPM   %Must cover the maximum operating flow for any pipe
PipesQLinMidPoints=[];
PipesHLinMidPoints=[];
PipesLinSegmSlope=[];
PipesLinSegmIntcpt=[];
for pp = 1 : PipeCount
    PLength=PipeLinkLength(pp);
    PDiameter=PipeLinkDiameter(pp)/ConfigurationConstants.FT2Inch;
    PRoughnessCoeff=PipeLinkRoughnessCoeff(pp); %HW Eq.
    [QLinMidPoints,HLinMidPoints,LinSegmSlope,LinSegmIntcpt]=PipeHeadLoss_Linearization(PLength,PDiameter,PRoughnessCoeff,Npw,QLinMax,ConfigurationConstants.GPMperCFS);
    PipesQLinMidPoints=[PipesQLinMidPoints, QLinMidPoints];
    PipesHLinMidPoints=[PipesHLinMidPoints, HLinMidPoints];
    PipesLinSegmSlope=[PipesLinSegmSlope, LinSegmSlope];
    PipesLinSegmIntcpt=[PipesLinSegmIntcpt, LinSegmIntcpt];
end

% Pump curve approximation
% Note that for networks with multiple Pumps you need to make this loop more
% dynamic and change the pump characteristics parameters
PumpCoeffs=[];
for MM=1: PumpCount
    PumpCurveApproximate;
    PumpCoeffs(:,MM)=PumpCoeff;
end
% Pump power consumption function relaxation and ensuring convexity
PumpPowerCoeffs=[];
for MM=1: PumpCount
    PumpPowerConsumptionApprox;
    PumpPowerCoeffs(:,MM)=PowerCoeff;
end

%% Initialize State Tracking Structures
TrueState = struct();
TrueState.TankHead = zeros(25, TankCount);
TrueState.JunctionHead = zeros(25, JunctionCount);
TrueState.PumpFlow = zeros(25, PumpCount);
TrueState.PumpSpeed = zeros(25, PumpCount);
TrueState.PipeFlow = zeros(25, PipeCount);
TrueState.PumpHead = zeros(25, PumpCount);
TrueState.ObjectiveValue = zeros(25, 1);

% Initial conditions (t=0)
TrueState.TankHead(1,:) = TankInitialHead;
TrueState.JunctionHead(1,:) = JunctionElevation;
TrueState.PipeFlow(1,:) = zeros(1, PipeCount);
TrueState.PumpFlow(1,:) = zeros(1, PumpCount);
TrueState.PumpSpeed(1,:) = zeros(1, PumpCount);
TrueState.PumpHead(1,:) = zeros(1, PumpCount);
TrueState.ObjectiveValue(1) = 0;

% Initialize Perceived State
PerceivedState = struct();
PerceivedState.TankHead = zeros(25, TankCount);
PerceivedState.JunctionHead = zeros(25, JunctionCount);
PerceivedState.PumpFlow = zeros(25, PumpCount);
PerceivedState.PumpSpeed = zeros(25, PumpCount);
PerceivedState.PipeFlow = zeros(25, PipeCount);
PerceivedState.PumpHead = zeros(25, PumpCount);
PerceivedState.ObjectiveValue = zeros(25, 1);

% Copy initial conditions to perceived state
PerceivedState.TankHead(1,:) = TrueState.TankHead(1,:);
PerceivedState.JunctionHead(1,:) = TrueState.JunctionHead(1,:);
PerceivedState.PipeFlow(1,:) = TrueState.PipeFlow(1,:);
PerceivedState.PumpFlow(1,:) = TrueState.PumpFlow(1,:);
PerceivedState.PumpSpeed(1,:) = TrueState.PumpSpeed(1,:);
PerceivedState.PumpHead(1,:) = TrueState.PumpHead(1,:);
PerceivedState.ObjectiveValue(1) = TrueState.ObjectiveValue(1);

% Attack parameters
targets.flow = 1;
targets.demand = 1;      % Junction 1 demand (only junction one has demand values in net1)
attack_magnitude = 0.1;  % Maximum 10% deviation

% Initialize tracking variables
EstimatedStatesOverTime = [];
MeasuredDataOverTime = [];
residualsOverTime = [];
FDI_data = [];
FDI_data_demand = [];
FDI_pump_data = [];
metric_data = [];
Final_data = [];

% ID parameters
n = 5;  % Number of measurements
CUSUM = zeros(n,25);

% The following CUSUM parameters are tuned for a desired false alarm rate
% of 1 false alarm in 24 hours using data during normal operation
mean_r = [0.5e-01  3.7722e-04  6.5746e-04   3.6390e-04   1.4357e-02]';
std_r =  [3.0033e-01  2.7758e-04  4.5796e-04   3.6549e-04   4.9949e-05]';
threshold = mean_r + 3*std_r;
bias = mean_r + 0.1*std_r;

% Define attack start and end time
att_start = 11;
att_end = 20;

% Define kwh price
kwh_price = 0.175;

%% Main Control Loop
tic
for tH = 1:24

    fprintf('Time step: %d\n', tH);

    %% True System Optimization
    % Define optimization variables
    JunctionsHeadVariable = sdpvar(1,JunctionCount);
    TanksHeadVariable = sdpvar(1,TankCount);
    PipesFlowVariable = sdpvar(1,PipeCount);
    PipesZetaVariable = sdpvar(Npw,PipeCount);
    PipesOmegaVariable = binvar(Npw,PipeCount);
    PumpsFlowVariable = sdpvar(1,PumpCount);
    PumpsHeadVariable = sdpvar(1,PumpCount);
    PumpSpeedVariable = sdpvar(1,PumpCount);

    HydOptObjectiveFunction = (PumpsFlowVariable/ConfigurationConstants.GPMperCFS)*PumpsHeadVariable*SG/(3960*eta_Pump)*0.746;

    % Basic constraints
    HydOptContraints = [];
    HydOptContraints = [HydOptContraints, JunctionsHeadVariable >= JunctionElevation];
    HydOptContraints = [HydOptContraints, PumpsFlowVariable >= zeros(PumpCount)];
    HydOptContraints = [HydOptContraints, PumpsHeadVariable >= zeros(PumpCount)];
    if TankCount > 0
        HydOptContraints = [HydOptContraints, TanksHeadVariable >= TankMinHead, TanksHeadVariable <= TankMaxHead];
    end

    % Mass balance constraints with true demands
    for JJ = 1:JunctionCount
        JuncConstTemp = JunctionDemand24(tH,JJ);
        for pp = 1:PipeCount
            if LinkFromTo(pp,2) == JJ
                JuncConstTemp = JuncConstTemp - PipesFlowVariable(1,pp);
            elseif LinkFromTo(pp,1) == JJ
                JuncConstTemp = JuncConstTemp + PipesFlowVariable(1,pp);
            end
        end
        if LinkFromTo(PumpIndex,2) == JJ
            JuncConstTemp = JuncConstTemp - PumpsFlowVariable;
            HydOptContraints = [HydOptContraints, ReservoirHead+PumpsHeadVariable == JunctionsHeadVariable(1,JJ)];
        elseif LinkFromTo(PumpIndex,1) == JJ
            JuncConstTemp = JuncConstTemp + PumpsFlowVariable;
            HydOptContraints = [HydOptContraints, -ReservoirHead+PumpsHeadVariable == -JunctionsHeadVariable(1,JJ)];
        end
        HydOptContraints = [HydOptContraints, JuncConstTemp == 0];
    end

    % Pipe dynamics for true system
    for pp = 1:PipeCount
        HydOptContraints = [HydOptContraints, PipesFlowVariable(1,pp)-ones(1,Npw)*PipesZetaVariable(:,pp) == 0];
        HydOptContraints = [HydOptContraints, ones(1,Npw)*PipesOmegaVariable(:,pp) == 1];

        for nn = 1:Npw
            HydOptContraints = [HydOptContraints, -PipesZetaVariable(nn,pp)+PipesQLinMidPoints(nn,pp)*PipesOmegaVariable(nn,pp) <= 0];
            HydOptContraints = [HydOptContraints, PipesZetaVariable(nn,pp)-PipesQLinMidPoints(nn+1,pp)*PipesOmegaVariable(nn,pp) <= 0];
        end

        if LinkFromTo(pp,1) <= JunctionCount && LinkFromTo(pp,2) <= JunctionCount
            JJn1 = LinkFromTo(pp,1);
            JJn2 = LinkFromTo(pp,2);
            HydOptContraints = [HydOptContraints, JunctionsHeadVariable(1,JJn1)-JunctionsHeadVariable(1,JJn2)-...
                PipesLinSegmSlope(:,pp)'*PipesZetaVariable(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable(:,pp) == 0];
        end

        if LinkFromTo (pp,1) == TankIndex
            HydOptContraints=[HydOptContraints, TanksHeadVariable==TrueState.TankHead(tH,1)-60*60/TankArea*(PipesFlowVariable(1,pp)/ConfigurationConstants.GPMperCFS)];
            JJn=LinkFromTo (pp,2);
            HydOptContraints=[HydOptContraints, -JunctionsHeadVariable(1,JJn)+TrueState.TankHead(tH,1)-PipesLinSegmSlope(:,pp)'*PipesZetaVariable(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable(:,pp)==0];
        elseif LinkFromTo (pp,2) == TankIndex
            HydOptContraints=[HydOptContraints, TanksHeadVariable==TrueState.TankHead(tH,1)+60*60/TankArea*(PipesFlowVariable(1,pp)/ConfigurationConstants.GPMperCFS)];
            JJn=LinkFromTo (pp,1);
            HydOptContraints=[HydOptContraints, JunctionsHeadVariable(1,JJn)-TrueState.TankHead(tH,1)-PipesLinSegmSlope(:,pp)'*PipesZetaVariable(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable(:,pp)==0];
        end
    end

    % Pump dynamics
    HydOptContraints = [HydOptContraints, PumpSpeedVariable >= 0, PumpSpeedVariable <= 1];
    HydOptContraints = [HydOptContraints, PumpsHeadVariable == PumpCoeff(1)*PumpsFlowVariable^2 + ...
        PumpCoeff(2)*PumpsFlowVariable + PumpCoeff(3)*PumpSpeedVariable^2 + ...
        PumpCoeff(4)*PumpSpeedVariable + PumpCoeff(5)];

    % Obj function for power is always positive
    HydOptContraints = [HydOptContraints, HydOptObjectiveFunction >= 0];

    % Solve true system optimization
    optimize(HydOptContraints, HydOptObjectiveFunction, sdpsettings('solver', 'gurobi', 'verbose', 0));

    % Store true state results
    if TankCount > 0
        TrueState.TankHead(tH+1,:)=value(TanksHeadVariable);
    end
    TrueState.JunctionHead(tH+1,:) = value(JunctionsHeadVariable);
    TrueState.PipeFlow(tH+1,:) = value(PipesFlowVariable);
    TrueState.PumpFlow(tH+1,:) = value(PumpsFlowVariable);
    TrueState.PumpSpeed(tH+1,:) = value(PumpSpeedVariable);
    TrueState.PumpHead(tH+1,:) = value(PumpsHeadVariable);
    TrueState.ObjectiveValue(tH+1) = value(HydOptObjectiveFunction)*ConfigurationConstants.GPMperCFS*kwh_price;

    %% Generate Measurements from normal operation
    flow_noise_std = 0.015;
    head_noise_std = 0.15;
    % Generate from true state
    [MeasuredFlows, MeasuredHeads] = generateNoisyMeasurements(TrueState.PipeFlow(tH+1,:), ...
        TrueState.JunctionHead(tH+1,:), flow_noise_std, head_noise_std);

    All_Data = [MeasuredFlows(:); MeasuredHeads(:)];
    num_flows = length(MeasuredFlows);
    num_heads = length(MeasuredHeads);

    if tH == 1
        total_measurements = length(All_Data);
        num_active_sensors = round((n/total_measurements) * total_measurements);

        % Chosen to ensure observability for Net1, 5 sensors (24% coverage)
        active_sensors =[1,5,7,9,15];

        % Create measurement matrix (fixed throughout simulation)
        MeasurementMatrix = zeros(num_active_sensors, total_measurements);
        for i = 1:num_active_sensors
            MeasurementMatrix(i, active_sensors(i)) = 1;
        end

        active_flow_sensors = active_sensors(active_sensors <= num_flows);
        active_head_sensors = active_sensors(active_sensors > num_flows);
        active_head_sensors_adjusted = active_head_sensors - num_flows;
    end

    MeasuredData = All_Data(active_sensors);

    %% FS-FDI Attack Implementation
    FS_FDI_active= 1;
    if FS_FDI_active == 1 && tH >= att_start && tH <= att_end
        % Create FDI variables
        FDI_vector = sdpvar(length(MeasuredData), 1);
        FDI_demand = sdpvar(1, 1);
        FDI_pump = sdpvar(1, 1);

        % Create attacked measurements
        attacked_measurements = MeasuredData + FDI_vector;
        attacked_flows = attacked_measurements(ismember(active_sensors, active_flow_sensors));

        % Attack constraints
        FDI_Constraints = [];

        % Target specific sensors
        flow_target_idx = find(active_sensors == targets.flow);
        for i = 1:length(MeasuredData)
            if i ~= flow_target_idx
                FDI_Constraints = [FDI_Constraints, FDI_vector(i) == 0];
            end
        end

        % Get current values
        current_demand = JunctionDemand24(tH, targets.demand);
        current_pipe_flow = TrueState.PipeFlow(tH+1, 1);  % Negative if flowing into junction
        current_pump_flow = TrueState.PumpFlow(tH+1);     % Positive into junction

        % Coordinated attack changes
        % 1. Increase demand
        FDI_Constraints = [FDI_Constraints, FDI_demand >= 0];
        FDI_Constraints = [FDI_Constraints, FDI_demand <= attack_magnitude * abs(current_demand)];

        % 2. Reduce pipe inflow (make less negative)
        flow_modification = FDI_vector(flow_target_idx);
        if current_pipe_flow < 0  % Flow into junction
            FDI_Constraints = [FDI_Constraints, flow_modification >= 0];  % Make flow less negative
            FDI_Constraints = [FDI_Constraints, flow_modification <= attack_magnitude * abs(current_pipe_flow)];
        end

        % 3. Pump modification to maintain mass balance
        % For mass balance: (pump_flow + FDI_pump) + (-(pipe_flow + flow_mod)) = demand + FDI_demand
        % Therefore: FDI_pump = FDI_demand - flow_modification
        FDI_Constraints = [FDI_Constraints, FDI_pump == FDI_demand - flow_modification];

        % Bound pump modification
        FDI_Constraints = [FDI_Constraints, abs(FDI_pump) <= attack_magnitude * abs(current_pump_flow)];

        % Mass balance verification
        if current_pipe_flow < 0
            mass_balance = (current_pump_flow + FDI_pump) + abs(current_pipe_flow + flow_modification) - (current_demand + FDI_demand);
        else
            mass_balance = (current_pump_flow + FDI_pump) - abs(current_pipe_flow + flow_modification) - (current_demand + FDI_demand);
        end
        FDI_Constraints = [FDI_Constraints, abs(mass_balance) <= 1e-6];

        H_meas = sqrt(W_meas)*H_direct;
        %-----------
        % 2) Build mass-balance constraints
        %-----------
        H_mass = zeros(JunctionCount, PipeCount + JunctionCount);
        z_mass = zeros(JunctionCount,1);

        for jJ = 1:JunctionCount
            dmd = JunctionDemand24(tH, jJ);
            if jj == targets.demand
                % If your attack changes the demand for that junction:
                z_mass(jJ) = dmd + value(FDI_demand);
            else
                z_mass(jJ) = dmd;
            end

            % sign convention: outflows are +, inflows are –
            for pp2 = 1:PipeCount
                if LinkFromTo(pp2,2) == jJ
                    H_mass(jJ, pp2) = -1;
                elseif LinkFromTo(pp2,1) == jJ
                    H_mass(jJ, pp2) = +1;
                end
            end
        end

        H_massW = H_mass;
        z_massW = z_mass;

        %-----------
        % 3) Build energy sub-block
        %-----------
        % a single slope from your piecewise set; for the pipe "pp2" => slopeVal
        % Then (HeadUp - HeadDown) - slopeVal*Flow = 0
        H_energy = [];
        z_energy = [];
        for pp2 = 1:PipeCount
            upNode = LinkFromTo(pp2,1);
            dnNode = LinkFromTo(pp2,2);
            if upNode <= JunctionCount && dnNode <= JunctionCount
                slopeVal = PipesLinSegmSlope(1, pp2); 
                rowVec = zeros(1, PipeCount + JunctionCount);
                rowVec(PipeCount+upNode) = +1;
                rowVec(PipeCount+dnNode) = -1;
                rowVec(pp2)            = -slopeVal;   
                H_energy = [H_energy; rowVec];
                z_energy = [z_energy; 0];
            end
        end
        H_energyW = H_energy;
        z_energyW = z_energy;

        H_system = [H_meas; H_massW; H_energyW];
        z_meas = sqrt(W_meas)*attacked_measurements;
        z = [z_meas; z_massW; z_energyW];

        EstimatedStates_FDI = inv(H_system.' * H_system)*H_system.'*z;

        attacked_residuals = attacked_measurements - MeasurementMatrix*EstimatedStates_FDI;
        FDI_Constraints = [FDI_Constraints, (CUSUM(flow_target_idx,tH) + abs(attacked_residuals(flow_target_idx)) - bias(flow_target_idx)) <= threshold(flow_target_idx)];

        % Maximize impact while maintaining feasibility
        optimize(FDI_Constraints, -(FDI_demand+FDI_vector(flow_target_idx)+FDI_pump), sdpsettings('solver', 'gurobi', 'verbose', 1));

        MeasuredData = MeasuredData + value(FDI_vector);
        FDI_data = [FDI_data, value(FDI_vector)];
        FDI_data_demand = [FDI_data_demand; value(FDI_demand)];
        FDI_pump_data = [FDI_pump_data; value(FDI_pump)];
    end

    %% State Estimation
    %-----------
    % 1) Build direct measurement sub-block
    %-----------
    numMeas = length(MeasuredData);
    H_direct = zeros(numMeas, PipeCount + JunctionCount);
    for iM = 1:numMeas
        globalIndex = active_sensors(iM);
        if globalIndex <= num_flows
            % Pipe flow sensor
            H_direct(iM, globalIndex) = 1;
        else
            % Head sensor => shift index by 'num_flows'
            juncIdx = globalIndex - num_flows;
            H_direct(iM, PipeCount + juncIdx) = 1;
        end
    end

    flow_weights = 1/(flow_noise_std^2)*ones(num_flows,1);
    head_weights = 1/(head_noise_std^2)*ones(num_heads,1);
    full_weights = [flow_weights; head_weights];
    activeWeights = full_weights(active_sensors);
    W_meas = diag(activeWeights);

    H_meas = sqrt(W_meas)*H_direct;           % Weighted measurement matrix
    z_meas = sqrt(W_meas)*MeasuredData;       % Weighted measurement vector

    %-----------
    % 2) Build mass-balance constraints
    %-----------
    H_mass = zeros(JunctionCount, PipeCount + JunctionCount);
    z_mass = zeros(JunctionCount,1);

    sigma_demand = 0.05;
    epsilon_demand = sigma_demand * randn(JunctionCount, 1);  % For each junction


    for jJ = 1:JunctionCount
        dmd = JunctionDemand24(tH, jJ);
        z_mass(jJ) = dmd;
        % sign convention: outflows are +, inflows are –
        for pp2 = 1:PipeCount
            if LinkFromTo(pp2,2) == jJ
                H_mass(jJ, pp2) = -1;
            elseif LinkFromTo(pp2,1) == jJ
                H_mass(jJ, pp2) = +1;
            end
        end
    end

    H_massW = H_mass;
    z_massW = z_mass;

    %-----------
    % 3) Build energy sub-block
    %-----------
    % a single slope from your piecewise set; for the pipe "pp2" => slopeVal
    % Then (HeadUp - HeadDown) - slopeVal*Flow = 0
    H_energy = [];
    z_energy = [];

    for pp2 = 1:PipeCount
        upNode = LinkFromTo(pp2,1);
        dnNode = LinkFromTo(pp2,2);
        if upNode <= JunctionCount && dnNode <= JunctionCount
            slopeVal = PipesLinSegmSlope(1, pp2);
            rowVec = zeros(1, PipeCount + JunctionCount);
            rowVec(PipeCount+upNode) = +1;
            rowVec(PipeCount+dnNode) = -1;
            rowVec(pp2)            = -slopeVal;     % minus slope * flow
            H_energy = [H_energy; rowVec];
            z_energy = [z_energy; 0];
        end
    end
    H_energyW = H_energy;
    z_energyW = z_energy;

    %-----------
    % 4) Combine everything
    %-----------
    H_system = [H_meas; H_massW; H_energyW];
    z_system = [z_meas; z_massW; z_energyW];

    %-----------
    % 5) Solve
    %-----------
    EstimatedStates = inv(H_system.' * H_system)*H_system.'*z_system;  % final estimate of [PipeFlows; JunctionHeads]

    % Store estimated states and measurements
    EstimatedStatesOverTime = [EstimatedStatesOverTime, EstimatedStates];
    MeasuredDataOverTime = [MeasuredDataOverTime, MeasuredData];

    % resid  = All_Data - EstimatedStates;
    % rmse   = sqrt(mean(resid.^2));
    % fprintf('Time step %d - WLS: final RMSE=%.4f\n', tH, rmse);


    %% Perceived System Optimization
    if tH >= att_start && tH <= att_end && FS_FDI_active == 1
        % Define optimization variables
        JunctionsHeadVariable_P = sdpvar(1,JunctionCount);
        TanksHeadVariable_P = sdpvar(1,TankCount);
        PipesFlowVariable_P = sdpvar(1,PipeCount);
        PipesZetaVariable_P = sdpvar(Npw,PipeCount);
        PipesOmegaVariable_P = binvar(Npw,PipeCount);
        PumpsFlowVariable_P = sdpvar(1,PumpCount);
        PumpsHeadVariable_P = sdpvar(1,PumpCount);
        PumpSpeedVariable_P = sdpvar(1,PumpCount);

        % Objective function
        HydOptObjectiveFunction_P = (PumpsFlowVariable_P/ConfigurationConstants.GPMperCFS)*PumpsHeadVariable_P*SG/(3960*eta_Pump)*0.746;

        % Constraints
        HydOptContraints_P = [];

        Input = MeasuredData';
        HydOptContraints_P = [HydOptContraints_P, abs(PipesFlowVariable_P(targets.flow)) <= abs(Input(flow_target_idx))];

        % Basic constraints
        HydOptContraints_P = [HydOptContraints_P, JunctionsHeadVariable_P >= JunctionElevation];
        HydOptContraints_P = [HydOptContraints_P, PumpsFlowVariable_P >= zeros(PumpCount)];
        HydOptContraints_P = [HydOptContraints_P, PumpsHeadVariable_P >= zeros(PumpCount)];
        if TankCount > 0
            HydOptContraints_P = [HydOptContraints_P, TanksHeadVariable_P >= TankMinHead, TanksHeadVariable_P <= TankMaxHead];
        end

        % Mass balance with attacked demands
        for JJ = 1:JunctionCount
            if JJ == targets.demand
                JuncConstTemp_P = JunctionDemand24(tH,JJ) + value(FDI_demand);
            else
                JuncConstTemp_P = JunctionDemand24(tH,JJ);
            end

            for pp = 1:PipeCount
                if LinkFromTo(pp,2) == JJ
                    JuncConstTemp_P = JuncConstTemp_P - PipesFlowVariable_P(1,pp);
                elseif LinkFromTo(pp,1) == JJ
                    JuncConstTemp_P = JuncConstTemp_P + PipesFlowVariable_P(1,pp);
                end
            end

            if LinkFromTo(PumpIndex,2) == JJ
                JuncConstTemp_P = JuncConstTemp_P - PumpsFlowVariable_P;
                HydOptContraints_P = [HydOptContraints_P, ReservoirHead+PumpsHeadVariable_P == JunctionsHeadVariable_P(1,JJ)];
            elseif LinkFromTo(PumpIndex,1) == JJ
                JuncConstTemp_P = JuncConstTemp_P + PumpsFlowVariable_P;
                HydOptContraints_P = [HydOptContraints_P, -ReservoirHead+PumpsHeadVariable_P == -JunctionsHeadVariable_P(1,JJ)];
            end

            HydOptContraints_P = [HydOptContraints_P, JuncConstTemp_P == 0];
        end

        % Pipe dynamics
        for pp = 1:PipeCount
            HydOptContraints_P = [HydOptContraints_P, PipesFlowVariable_P(1,pp)-ones(1,Npw)*PipesZetaVariable_P(:,pp) == 0];
            HydOptContraints_P = [HydOptContraints_P, ones(1,Npw)*PipesOmegaVariable_P(:,pp) == 1];

            for nn = 1:Npw
                HydOptContraints_P = [HydOptContraints_P, -PipesZetaVariable_P(nn,pp)+PipesQLinMidPoints(nn,pp)*PipesOmegaVariable_P(nn,pp) <= 0];
                HydOptContraints_P = [HydOptContraints_P, PipesZetaVariable_P(nn,pp)-PipesQLinMidPoints(nn+1,pp)*PipesOmegaVariable_P(nn,pp) <= 0];
            end

            if LinkFromTo(pp,1) <= JunctionCount && LinkFromTo(pp,2) <= JunctionCount
                JJn1 = LinkFromTo(pp,1);
                JJn2 = LinkFromTo(pp,2);
                HydOptContraints_P = [HydOptContraints_P, JunctionsHeadVariable_P(1,JJn1)-JunctionsHeadVariable_P(1,JJn2)-...
                    PipesLinSegmSlope(:,pp)'*PipesZetaVariable_P(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable_P(:,pp) == 0];
            end

            % Tank connections
            if LinkFromTo(pp,1) == TankIndex
                HydOptContraints_P = [HydOptContraints_P, TanksHeadVariable_P==PerceivedState.TankHead(tH,1)-60*60/TankArea*(PipesFlowVariable_P(1,pp)/ConfigurationConstants.GPMperCFS)];
                JJn = LinkFromTo(pp,2);
                HydOptContraints_P = [HydOptContraints_P, -JunctionsHeadVariable_P(1,JJn)+PerceivedState.TankHead(tH,1)-PipesLinSegmSlope(:,pp)'*PipesZetaVariable_P(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable_P(:,pp)==0];
            elseif LinkFromTo(pp,2) == TankIndex
                HydOptContraints_P = [HydOptContraints_P, TanksHeadVariable_P==PerceivedState.TankHead(tH,1)+60*60/TankArea*(PipesFlowVariable_P(1,pp)/ConfigurationConstants.GPMperCFS)];
                JJn = LinkFromTo(pp,1);
                HydOptContraints_P = [HydOptContraints_P, JunctionsHeadVariable_P(1,JJn)-PerceivedState.TankHead(tH,1)-PipesLinSegmSlope(:,pp)'*PipesZetaVariable_P(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable_P(:,pp)==0];
            end
        end

        % Pump dynamics
        HydOptContraints_P = [HydOptContraints_P, PumpSpeedVariable_P >= 0, PumpSpeedVariable_P <= 1];
        HydOptContraints_P = [HydOptContraints_P, PumpsHeadVariable_P == PumpCoeff(1)*PumpsFlowVariable_P^2 + ...
            PumpCoeff(2)*PumpsFlowVariable_P + PumpCoeff(3)*PumpSpeedVariable_P^2 + ...
            PumpCoeff(4)*PumpSpeedVariable_P + PumpCoeff(5)];

        % Obj function for power is always positive
        HydOptContraints_P = [HydOptContraints_P, HydOptObjectiveFunction_P >= 0];

    else
        % During normal operation - same as true system
        JunctionsHeadVariable_P = sdpvar(1,JunctionCount);
        TanksHeadVariable_P = sdpvar(1,TankCount);
        PipesFlowVariable_P = sdpvar(1,PipeCount);
        PipesZetaVariable_P = sdpvar(Npw,PipeCount);
        PipesOmegaVariable_P = binvar(Npw,PipeCount);
        PumpsFlowVariable_P = sdpvar(1,PumpCount);
        PumpsHeadVariable_P = sdpvar(1,PumpCount);
        PumpSpeedVariable_P = sdpvar(1,PumpCount);

        HydOptObjectiveFunction_P = (PumpsFlowVariable_P/ConfigurationConstants.GPMperCFS)*PumpsHeadVariable_P*SG/(3960*eta_Pump)*0.746;

        HydOptContraints_P = [];
        HydOptContraints_P = [HydOptContraints_P, JunctionsHeadVariable_P >= JunctionElevation];
        HydOptContraints_P = [HydOptContraints_P, PumpsFlowVariable_P >= zeros(PumpCount)];
        HydOptContraints_P = [HydOptContraints_P, PumpsHeadVariable_P >= zeros(PumpCount)];
        if TankCount > 0
            HydOptContraints_P = [HydOptContraints_P, TanksHeadVariable_P >= TankMinHead, TanksHeadVariable_P <= TankMaxHead];
        end

        for JJ = 1:JunctionCount
            JuncConstTemp_P = JunctionDemand24(tH,JJ);
            for pp = 1:PipeCount
                if LinkFromTo(pp,2) == JJ
                    JuncConstTemp_P = JuncConstTemp_P - PipesFlowVariable_P(1,pp);
                elseif LinkFromTo(pp,1) == JJ
                    JuncConstTemp_P = JuncConstTemp_P + PipesFlowVariable_P(1,pp);
                end
            end

            if LinkFromTo(PumpIndex,2) == JJ
                JuncConstTemp_P = JuncConstTemp_P - PumpsFlowVariable_P;
                HydOptContraints_P = [HydOptContraints_P, ReservoirHead+PumpsHeadVariable_P == JunctionsHeadVariable_P(1,JJ)];
            elseif LinkFromTo(PumpIndex,1) == JJ
                JuncConstTemp_P = JuncConstTemp_P + PumpsFlowVariable_P;
                HydOptContraints_P = [HydOptContraints_P, -ReservoirHead+PumpsHeadVariable_P == -JunctionsHeadVariable_P(1,JJ)];
            end

            HydOptContraints_P = [HydOptContraints_P, JuncConstTemp_P == 0];
        end

        for pp = 1:PipeCount
            HydOptContraints_P = [HydOptContraints_P, PipesFlowVariable_P(1,pp)-ones(1,Npw)*PipesZetaVariable_P(:,pp) == 0];
            HydOptContraints_P = [HydOptContraints_P, ones(1,Npw)*PipesOmegaVariable_P(:,pp) == 1];

            for nn = 1:Npw
                HydOptContraints_P = [HydOptContraints_P, -PipesZetaVariable_P(nn,pp)+PipesQLinMidPoints(nn,pp)*PipesOmegaVariable_P(nn,pp) <= 0];
                HydOptContraints_P = [HydOptContraints_P, PipesZetaVariable_P(nn,pp)-PipesQLinMidPoints(nn+1,pp)*PipesOmegaVariable_P(nn,pp) <= 0];
            end

            if LinkFromTo(pp,1) <= JunctionCount && LinkFromTo(pp,2) <= JunctionCount
                JJn1 = LinkFromTo(pp,1);
                JJn2 = LinkFromTo(pp,2);
                HydOptContraints_P = [HydOptContraints_P, JunctionsHeadVariable_P(1,JJn1)-JunctionsHeadVariable_P(1,JJn2)-...
                    PipesLinSegmSlope(:,pp)'*PipesZetaVariable_P(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable_P(:,pp) == 0];
            end

            if LinkFromTo (pp,1) == TankIndex
                HydOptContraints_P=[HydOptContraints_P, TanksHeadVariable_P==PerceivedState.TankHead(tH,1)-60*60/TankArea*(PipesFlowVariable_P(1,pp)/ConfigurationConstants.GPMperCFS)];
                JJn=LinkFromTo (pp,2);
                HydOptContraints_P=[HydOptContraints_P, -JunctionsHeadVariable_P(1,JJn)+PerceivedState.TankHead(tH,1)-PipesLinSegmSlope(:,pp)'*PipesZetaVariable_P(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable_P(:,pp)==0];
            elseif LinkFromTo (pp,2) == TankIndex
                HydOptContraints_P=[HydOptContraints_P, TanksHeadVariable_P==PerceivedState.TankHead(tH,1)+60*60/TankArea*(PipesFlowVariable_P(1,pp)/ConfigurationConstants.GPMperCFS)];
                JJn=LinkFromTo (pp,1);
                HydOptContraints_P=[HydOptContraints_P, JunctionsHeadVariable_P(1,JJn)-PerceivedState.TankHead(tH,1)-PipesLinSegmSlope(:,pp)'*PipesZetaVariable_P(:,pp)-PipesLinSegmIntcpt(:,pp)'*PipesOmegaVariable_P(:,pp)==0];
            end
        end

        HydOptContraints_P = [HydOptContraints_P, PumpSpeedVariable_P >= 0, PumpSpeedVariable_P <= 1];
        HydOptContraints_P = [HydOptContraints_P, PumpsHeadVariable_P == PumpCoeff(1)*PumpsFlowVariable_P^2 + ...
            PumpCoeff(2)*PumpsFlowVariable_P + PumpCoeff(3)*PumpSpeedVariable_P^2 + ...
            PumpCoeff(4)*PumpSpeedVariable_P + PumpCoeff(5)];

        HydOptContraints_P = [HydOptContraints_P, HydOptObjectiveFunction_P >= 0];
    end

    % Solve perceived system optimization
    optimize(HydOptContraints_P, HydOptObjectiveFunction_P, sdpsettings('solver', 'gurobi', 'verbose', 0));

    % Store perceived system results
    if TankCount > 0
        PerceivedState.TankHead(tH+1,:) = value(TanksHeadVariable_P);
    end
    PerceivedState.JunctionHead(tH+1,:) = value(JunctionsHeadVariable_P);
    PerceivedState.PipeFlow(tH+1,:) = value(PipesFlowVariable_P);
    PerceivedState.PumpFlow(tH+1,:) = value(PumpsFlowVariable_P);
    PerceivedState.PumpSpeed(tH+1,:) = value(PumpSpeedVariable_P);
    PerceivedState.PumpHead(tH+1,:) = value(PumpsHeadVariable_P);
    PerceivedState.ObjectiveValue(tH+1) = value(HydOptObjectiveFunction_P)*ConfigurationConstants.GPMperCFS*kwh_price;

    %% Intrusion Detection
    residuals = MeasuredData - MeasurementMatrix*EstimatedStates;
    residualsOverTime = [residualsOverTime, residuals];
    flow_sensor_idx = find(active_sensors == targets.flow);
    CUSUM_metric = abs(residuals);
    CUSUM(:,tH+1) = max(CUSUM(:,tH) + CUSUM_metric - bias, 0);

    % Store alarm times and values
    if tH == 1
        alarm_times = [];
        alarm_values = [];
    end

    % Check for alarms
    for i = 1:n
        if CUSUM(i, tH+1) > threshold(i)
            fprintf('Alarm triggered by measurement %d at hour = %d, value = %f\n', ...
                i, tH, CUSUM(i, tH+1));
            % Store alarm info before reset
            if i == flow_sensor_idx
                alarm_times = [alarm_times, tH];
                alarm_values = [alarm_values, CUSUM(i, tH+1)];
            end
            CUSUM(i, tH+1) = 0;
        end
    end

    % Store metrics
    metric_data = [metric_data, abs(MeasuredData - MeasurementMatrix*EstimatedStates)];
    Final_data = [Final_data, MeasuredData];

end  % End of main time loop
Time_OptProb = toc;

%% Pseudo pump speed for true data
% Checking for psuedo pump speed
% For reasons why this might happen refer to the paper
% "Quality-Aware Hydraulic Control in Drinking Water Networks via Controllability Proxies"
HydOptOutObjFuncSeries=[0];
PumpStatus=cell(25,1);
PumpStatus{1,1}= '-';
for tH = 1 : 24
    if PerceivedState.PumpFlow(tH+1)<1 || PerceivedState.PumpHead(tH+1)<1
        PerceivedState.PumpSpeed(tH+1)=0;
        HydOptOutObjFuncSeries(tH+1)=0;
        PumpStatus{tH+1,1}= 'Off';
    else
        PumpStatus{tH+1,1}= 'On';
    end
end
%% Pseudo pump speed for precieved data
HydOptOutObjFuncSeries=[0];
PumpStatus=cell(25,1);
PumpStatus{1,1}= '-';
for tH = 1 : 24
    if TrueState.PumpFlow(tH+1)<1 || TrueState.PumpHead(tH+1)<1
        TrueState.PumpSpeed(tH+1)=0;
        HydOptOutObjFuncSeries(tH+1)=0;
        PumpStatus{tH+1,1}= 'Off';
    else
        PumpStatus{tH+1,1}= 'On';
    end
end

%% Save results for analysis
results = struct();
results.TrueState = TrueState;
results.PerceivedState = PerceivedState;
results.EstimatedStatesOverTime = EstimatedStatesOverTime;
results.MeasuredDataOverTime = MeasuredDataOverTime;
results.residualsOverTime = residualsOverTime;
results.FDI_data = FDI_data;
results.FDI_data_demand = FDI_data_demand;
results.CUSUM = CUSUM;
results.Time_OptProb = Time_OptProb;

save('attack_simulation_results.mat', 'results');

%% Plotting
plot_active = 1;
if plot_active == 1

    %% Clear workspace and set up plotting defaults
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'DefaultAxesFontSize', 11);

    % Load simulation results if they exist
    if exist('attack_simulation_results.mat', 'file')
        load('attack_simulation_results.mat');
    end

    %% Plot 1: Mass Balance Validation
    flow_target_idx = find(active_sensors == targets.flow);
    FDI_flow = FDI_data(flow_target_idx,:)';
    FDI_flow_all = zeros(24,1);
    FDI_flow_all(att_start:att_end) = FDI_flow;

    mass_errors = zeros(1, 24);
    for t = 1:24
        fprintf('\nTime step %d:\n', t);

        % Get base flows
        base_pipe_flow = TrueState.PipeFlow(t+1, 1); % Keep original sign
        base_pump_flow = TrueState.PumpFlow(t+1); % Always into junction

        % Get modified flows during attack
        pipe_flow = base_pipe_flow;
        pump_flow = base_pump_flow;
        if t >= att_start && t <= att_end
            attack_idx = t - att_start + 1;
            pipe_flow = base_pipe_flow + FDI_flow_all(t);
            pump_flow = base_pump_flow + FDI_pump_data(attack_idx);
        end

        % Get demand
        demand = JunctionDemand24(t,1);
        if t >= att_start && t <= att_end
            attack_idx = t - att_start + 1;
            demand = demand + FDI_data_demand(attack_idx);
        end

        % Calculate mass balance based on flow direction
        if pipe_flow < 0 % Pipe flow into junction
            total_inflow = pump_flow + abs(pipe_flow);
        else % Pipe flow out of junction
            total_inflow = pump_flow - abs(pipe_flow);
        end

        fprintf('Base pump flow: %.6f\n', base_pump_flow);
        fprintf('Modified pump flow: %.6f\n', pump_flow);
        fprintf('Base pipe flow: %.6f\n', base_pipe_flow);
        fprintf('Modified pipe flow: %.6f\n', pipe_flow);
        fprintf('Total inflow: %.6f\n', total_inflow);
        fprintf('Demand: %.6f\n', demand);

        % Mass balance error
        mass_errors(t) = abs(total_inflow - demand);
        fprintf('Mass balance error: %.6f\n', mass_errors(t));
    end

    % Plot validation
    figure('Position', [100 100 800 400]);
    semilogy(1:24, mass_errors, 'k-', 'LineWidth', 2);
    hold on;
    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    xlabel('Time (hours)');
    ylabel('Mass Balance Error (GPM)');
    title('\textbf{Mass Balance Validation}');
    grid on;

    %% Plot 2: Hydraulic Impact Analysis
    figure('Position', [100 100 1200 400]);

    % Create new figure and close any leftover legend handles
    clf;
    legend('off');

    % Create proper time vector
    t = 1:24;  % Since we're plotting hourly data

    % 1. Pump Speed Response
    subplot(1,4,1)
    hold on
    plot(t, TrueState.PumpSpeed(2:end), 'b-', 'LineWidth', 2, 'DisplayName', 'True Speed');
    plot(t, PerceivedState.PumpSpeed(2:end), 'r--', 'LineWidth', 2, 'DisplayName', 'Under Attack');
    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(att_start, 'k--', 'HandleVisibility', 'off');
    xline(att_end, 'k--', 'HandleVisibility', 'off');
    xlabel('Time (hours)');
    ylabel('Relative Speed');
    title('\textbf{(a) Pump Speed Response}');
    legend('Location', 'best');
    grid on;

    zoom_start = 10;
    zoom_end = 14;

    % 2. Pump Flow Response
    subplot(1,4,2)
    hold on
    plot(t, TrueState.PumpFlow(2:end), 'b-', 'LineWidth', 2, 'DisplayName', 'True Flow');
    plot(t, PerceivedState.PumpFlow(2:end), 'r--', 'LineWidth', 2, 'DisplayName', 'Under Attack');
    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(att_start, 'k--', 'HandleVisibility', 'off');
    xline(att_end, 'k--', 'HandleVisibility', 'off');
    xlabel('Time (hours)');
    ylabel('Flow Rate (GPM)');
    title('\textbf{(b) Pump Flow Response}');
    legend('Location', 'northwest');
    grid on;

    % 2. Pump Head
    subplot(1,4,3)
    hold on
    plot(t, TrueState.PumpHead(2:end), 'b-', 'LineWidth', 2, 'DisplayName', 'True Head');
    plot(t, PerceivedState.PumpHead(2:end), 'r--', 'LineWidth', 2, 'DisplayName', 'Under Attack');
    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(att_start, 'k--', 'HandleVisibility', 'off');
    xline(att_end, 'k--', 'HandleVisibility', 'off');
    xlabel('Time (hours)');
    ylabel('Head (ft)');
    title('\textbf{(b) Pump Head}');
    legend('Location', 'best');
    grid on;

    % 3. Hourly Operating Cost
    subplot(1,4,4)
    hold on
    % Prepare data for bar plot
    cost_data = zeros(24,2);
    cost_data(:,1) = TrueState.ObjectiveValue(2:end);
    cost_data(:,2) = PerceivedState.ObjectiveValue(2:end);

    % Create grouped bar plot
    b = bar(t, cost_data, 'grouped');
    % b(1).FaceColor = [0 0 1]; % Blue for true cost
    % b(2).FaceColor = [1 0 0]; % Red for perceived cost

    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(att_start, 'k--', 'HandleVisibility', 'off');
    xline(att_end, 'k--', 'HandleVisibility', 'off');
    xlabel('Time (hours)');
    ylabel('Cost (\$/hour)');
    title('\textbf{(c) Hourly Operating Cost}');
    legend('True Cost', 'Under Attack', 'Location', 'best');
    grid on;

    % sgtitle('\textbf{Hydraulic Impact Analysis}', 'FontSize', 14);
    set(gcf, 'Color', 'white');

    % Save high-quality figure
    exportgraphics(gcf, 'hydraulic_impact.pdf', 'ContentType', 'vector');

    %% Plot 3: Pump flow, attack values, and CUSUM statistic
    % Replace Pipe Flow with Pump Flow in the main plot
    t=0:24;
    subplot(3,1,1)
    hold on
    plot(t, TrueState.PumpFlow(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'True Pump Flow');
    plot(t, PerceivedState.PumpFlow(:,1), 'r--', 'LineWidth', 2, 'DisplayName', 'Under Attack');
    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(att_start, 'k--', 'HandleVisibility', 'off');
    xline(att_end, 'k--', 'HandleVisibility', 'off');
    xlabel('Time (hours)');
    ylabel('Flow Rate (GPM)');
    title('\textbf{(b) Pump Flow Response}');
    legend('Location', 'northwest');
    grid on;

    % Add a new plot to show attack on demand, pipe flow, and pump flow
    subplot(3,1,2)
    attack_hours = att_start:att_end;
    bar_data = [FDI_data(flow_target_idx,:)', FDI_data_demand, FDI_pump_data];
    b = bar(attack_hours, bar_data, 'stacked');
    b(1).FaceColor = [0.2 0.2 0.8]; % Pipe Flow Attack
    b(2).FaceColor = [0.8 0.2 0.2]; % Demand Attack
    b(3).FaceColor = [0.2 0.8 0.2]; % Pump Flow Attack
    xlabel('Time (hours)');
    ylabel('Attack Magnitude');
    title('\textbf{(d) Attack Components: Pipe, Demand, and Pump}');
    legend('Pipe Flow Attack', 'Demand Attack', 'Pump Flow Attack', 'Location', 'northwest');
    grid on;

    % Filter alarms to only include those within the attack period
    valid_alarm_indices = alarm_times >= att_start & alarm_times <= att_end;
    filtered_alarm_times = alarm_times(valid_alarm_indices);
    filtered_alarm_values = alarm_values(valid_alarm_indices);

    % CUSUM Detection Statistics (unchanged for context)
    subplot(3,1,3)
    hold on
    flow_sensor_idx = find(active_sensors == targets.flow);
    plot(0:24, CUSUM(flow_sensor_idx,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CUSUM Statistic');
    plot(0:24, threshold(flow_sensor_idx)*ones(1,25), 'r--', 'LineWidth', 1.5, 'DisplayName', 'CUSUM Threshold');
    if ~isempty(filtered_alarm_times)
        scatter(filtered_alarm_times, filtered_alarm_values, 100, 'r', 'filled', 'diamond', ...
            'DisplayName', 'Alarms', 'HandleVisibility', 'on');
        for i = 1:length(filtered_alarm_times)
            arrow_x = filtered_alarm_times(i);
            arrow_y = [filtered_alarm_values(i), 0];
            plot([arrow_x, arrow_x], arrow_y, 'r:', 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(att_start, 'k--', 'HandleVisibility', 'off');
    xline(att_end, 'k--', 'HandleVisibility', 'off');
    xlabel('Time (hours)');
    ylabel('CUSUM Statistic');
    title('\textbf{(e) CUSUM Detection Statistics}');
    legend('Location', 'northwest');
    grid on;

    %% Plot 4: Operating Points and Cost Analysis
    figure('Position', [100 100 1200 400]);

    % 1. Operating Point Trajectory - Clear Version
    subplot(2,1,1)
    hold on

    % Plot points for normal operation (before attack)
    scatter(TrueState.PumpFlow(2:att_start), TrueState.PumpHead(2:att_start), 50, 'b', 'filled', 'DisplayName', 'Normal Operation');

    % Plot points during attack period
    scatter(TrueState.PumpFlow(att_start+1:att_end+1), TrueState.PumpHead(att_start+1:att_end+1), 50, 'b^', 'filled', 'DisplayName', 'True (During Attack)');
    scatter(PerceivedState.PumpFlow(att_start+1:att_end+1), PerceivedState.PumpHead(att_start+1:att_end+1), 50, 'r^', 'filled', 'DisplayName', 'Under Attack');

    % Plot points after attack
    scatter(TrueState.PumpFlow(att_end+2:end), TrueState.PumpHead(att_end+2:end), 50, 'b', 'filled', 'HandleVisibility', 'off');

    % Add pump curves for reference
    Q = 0:100:1000;
    speeds = [0.6 0.8 1.0]; % Example pump speeds
    for s = speeds
        H = s^2 * (PumpCoeff(3)) + s*PumpCoeff(4) + PumpCoeff(5) + ...
            PumpCoeff(1)*Q.^2 + PumpCoeff(2)*Q;
        plot(Q, H, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
        % Add speed label
        text(Q(end), H(end), sprintf('%.1fx', s), 'FontSize', 8);
    end

    xlabel('Flow Rate (GPM)');
    ylabel('Head (ft)');
    title('\textbf{(a) Pump Operating Points}');
    legend('Location', 'best');
    grid on;

    set(gcf, 'Color', 'white');
    % 2. Cumulative Cost Impact
    subplot(2,1,2)
    hold on

    % Calculate baseline cumulative cost
    cumul_cost_true = cumsum(TrueState.ObjectiveValue);
    cumul_cost_attack = cumsum(PerceivedState.ObjectiveValue);

    % Calculate percentage increase starting from first significant cost (hour 15)
    cost_increase_percent = zeros(size(cumul_cost_true));
    start_hour = 15;  % When we see first significant cost
    cost_increase_percent(start_hour+1:end) = ...
        (cumul_cost_attack(start_hour+1:end) - cumul_cost_true(start_hour+1:end)) ./ ...
        cumul_cost_true(start_hour+1:end) * 100;

    % Plot absolute costs on left y-axis
    yyaxis left
    ax = gca;
    ax.YColor = 'k';
    plot(0:24, cumul_cost_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Cost');
    plot(0:24, cumul_cost_attack, 'r--', 'LineWidth', 2, 'DisplayName', 'Under Attack');
    ylabel('Cumulative Cost (\$)');

    % Plot percentage increase on right y-axis
    yyaxis right
    ax.YColor = [0.2 0.2 0.2];
    plot(0:24, cost_increase_percent, 'Color', [0.2 0.2 0.2], 'LineStyle', ':', ...
        'LineWidth', 2, 'DisplayName', 'Increase (\%)');
    ylabel('Cost Increase (\%)');
    ylim([0 100]);  % Limit y-axis to reasonable percentage range

    % Add attack period shading
    y_limits = ylim;
    patch([att_start att_end att_end att_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(att_start, 'k--', 'HandleVisibility', 'off');
    xline(att_end, 'k--', 'HandleVisibility', 'off');

    xlabel('Time (hours)');
    title('\textbf{(b) Cumulative Operating Cost Impact}');
    legend('Location', 'northwest');
    grid on;
end

%% Generating noisy measurements function
function [MeasuredFlows, MeasuredHeads] = generateNoisyMeasurements(TruePipeFlow, TrueJunctionHead, flow_noise_std, head_noise_std)
% Set the random seed for reproducible noise
rng(1);  % Change the seed to get different fixed noise patterns

% Generate fixed noise for flows and heads
flow_noise = flow_noise_std * randn(size(TruePipeFlow)) + 0.4 * rand(size(TruePipeFlow));
head_noise = head_noise_std * randn(size(TrueJunctionHead)) + 0.4 * rand(size(TrueJunctionHead));

% Generate noisy measurements based on the provided noise standard deviations
flow_noise(1)=0;
MeasuredFlows = TruePipeFlow + flow_noise;
MeasuredHeads = TrueJunctionHead + head_noise;
end