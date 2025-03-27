%% Network Analysis Script
function analyzeNetwork(G)
    % Get basic network information
    node_id = G.getNodeNameID;
    PipeIndex = G.getLinkPipeIndex;
    PumpIndex = G.getLinkPumpIndex;
    JunctionIndex = G.getNodeJunctionIndex;
    ReservoirIndex = G.getNodeReservoirIndex;
    TankIndex = G.getNodeTankIndex;
    LinkFromTo = G.getLinkNodesIndex;
    
    fprintf('\n=== Network Overview ===\n');
    fprintf('Junctions: %d\n', length(JunctionIndex));
    fprintf('Pipes: %d\n', length(PipeIndex));
    fprintf('Pumps: %d\n', length(PumpIndex));
    fprintf('Tanks: %d\n', length(TankIndex));
    fprintf('Reservoirs: %d\n', length(ReservoirIndex));
    
    % Analyze pump connections
    fprintf('\n=== Pump Information ===\n');
    for p = 1:length(PumpIndex)
        from_node = LinkFromTo(PumpIndex(p), 1);
        to_node = LinkFromTo(PumpIndex(p), 2);
        fprintf('Pump %d: From node %d to node %d\n', p, from_node, to_node);
        
        % Find pipes connected to these nodes
        connected_pipes_from = find(LinkFromTo(:,1) == from_node | LinkFromTo(:,2) == from_node);
        connected_pipes_to = find(LinkFromTo(:,1) == to_node | LinkFromTo(:,2) == to_node);
        
        fprintf('  Connected pipes at source: %s\n', mat2str(connected_pipes_from));
        fprintf('  Connected pipes at destination: %s\n', mat2str(connected_pipes_to));
    end
    
    % Identify critical junctions (nodes with multiple connections)
    fprintf('\n=== Critical Junctions ===\n');
    for j = 1:length(JunctionIndex)
        node = JunctionIndex(j);
        connected_pipes = find(LinkFromTo(:,1) == node | LinkFromTo(:,2) == node);
        if length(connected_pipes) > 2  % More than 2 connections
            fprintf('Junction %d has %d connections. Connected pipes: %s\n', ...
                node, length(connected_pipes), mat2str(connected_pipes));
            
            % Check if connected to any pumps
            connected_pumps = find(LinkFromTo(PumpIndex,1) == node | LinkFromTo(PumpIndex,2) == node);
            if ~isempty(connected_pumps)
                fprintf('  Connected to pump(s): %s\n', mat2str(connected_pumps));
            end
        end
    end
    
    % Find potential attack targets
    fprintf('\n=== Potential Attack Targets ===\n');
    for p = 1:length(PipeIndex)
        node1 = LinkFromTo(p, 1);
        node2 = LinkFromTo(p, 2);
        
        % Count connections at each end
        pipes_node1 = find(LinkFromTo(:,1) == node1 | LinkFromTo(:,2) == node1);
        pipes_node2 = find(LinkFromTo(:,1) == node2 | LinkFromTo(:,2) == node2);
        
        % Check if near pumps
        near_pump1 = any(LinkFromTo(PumpIndex,1) == node1 | LinkFromTo(PumpIndex,2) == node1);
        near_pump2 = any(LinkFromTo(PumpIndex,1) == node2 | LinkFromTo(PumpIndex,2) == node2);
        
        % If pipe connects to a junction with multiple connections and near a pump
        if (length(pipes_node1) >= 2 || length(pipes_node2) >= 2) && (near_pump1 || near_pump2)
            fprintf('Potential target - Pipe %d:\n', p);
            fprintf('  Connects nodes %d (%d connections) and %d (%d connections)\n', ...
                node1, length(pipes_node1), node2, length(pipes_node2));
            if near_pump1
                pump_idx = find(LinkFromTo(PumpIndex,1) == node1 | LinkFromTo(PumpIndex,2) == node1);
                fprintf('  Node %d connected to pump(s): %s\n', node1, mat2str(pump_idx));
            end
            if near_pump2
                pump_idx = find(LinkFromTo(PumpIndex,1) == node2 | LinkFromTo(PumpIndex,2) == node2);
                fprintf('  Node %d connected to pump(s): %s\n', node2, mat2str(pump_idx));
            end
        end
    end
    
    % Analyze base demands
    fprintf('\n=== Demand Analysis ===\n');
    BaseDemand = G.getNodeBaseDemands{1};
    BaseDemand = BaseDemand(JunctionIndex);
    
    [sorted_demands, demand_idx] = sort(BaseDemand, 'descend');
    fprintf('Top 10 demand nodes:\n');
    for i = 1:min(10, length(sorted_demands))
        node = JunctionIndex(demand_idx(i));
        fprintf('Node %d: %.2f\n', node, sorted_demands(i));
    end
end