function checkNetworkBalance(TrueState, JunctionDemand24, LinkFromTo, JunctionCount, PipeCount, PumpCount, PumpIndex)
    % Initialize storage for errors
    mass_balance_errors = zeros(24, JunctionCount);
    energy_balance_errors = zeros(24, PipeCount);
    
    fprintf('Checking network balance for 24 hours...\n');
    
    for t = 1:24
        % Check mass balance at each junction
        for j = 1:JunctionCount
            % Initialize balance for this junction
            balance = -JunctionDemand24(t,j);  % Start with negative demand
            
            % Get all pipes connected to this junction
            pipe_indices_to = find(LinkFromTo(:,2) == j);  % Flows into junction
            pipe_indices_from = find(LinkFromTo(:,1) == j); % Flows out of junction
            
            % Add inflows
            for idx = pipe_indices_to'
                if idx <= PipeCount  % Check if it's a pipe (not a pump)
                    balance = balance + TrueState.PipeFlow(t+1,idx);
                end
            end
            
            % Subtract outflows
            for idx = pipe_indices_from'
                if idx <= PipeCount  % Check if it's a pipe (not a pump)
                    balance = balance - TrueState.PipeFlow(t+1,idx);
                end
            end
            
            % Add pump contributions
            for p = 1:PumpCount
                pump_idx = PumpIndex(p);
                if LinkFromTo(pump_idx,2) == j  % Pump discharge into junction
                    balance = balance + TrueState.PumpFlow(t+1,p);
                elseif LinkFromTo(pump_idx,1) == j  % Pump suction from junction
                    balance = balance - TrueState.PumpFlow(t+1,p);
                end
            end
            
            % Store absolute error
            mass_balance_errors(t,j) = abs(balance);
            
            % Debug output for large imbalances
            if abs(balance) > 1.0  % Threshold of 1 GPM
                fprintf('\nLarge imbalance at junction %d, time %d: %.6f GPM\n', j, t, balance);
                fprintf('  Demand: %.6f GPM\n', JunctionDemand24(t,j));
                fprintf('  Connected pipes in: %s\n', mat2str(pipe_indices_to));
                fprintf('  Connected pipes out: %s\n', mat2str(pipe_indices_from));
                for idx = pipe_indices_to'
                    if idx <= PipeCount
                        fprintf('  Inflow pipe %d: %.6f GPM\n', idx, TrueState.PipeFlow(t+1,idx));
                    end
                end
                for idx = pipe_indices_from'
                    if idx <= PipeCount
                        fprintf('  Outflow pipe %d: %.6f GPM\n', idx, TrueState.PipeFlow(t+1,idx));
                    end
                end
                % Check for connected pumps
                for p = 1:PumpCount
                    pump_idx = PumpIndex(p);
                    if LinkFromTo(pump_idx,2) == j
                        fprintf('  Pump %d inflow: %.6f GPM\n', p, TrueState.PumpFlow(t+1,p));
                    elseif LinkFromTo(pump_idx,1) == j
                        fprintf('  Pump %d outflow: %.6f GPM\n', p, TrueState.PumpFlow(t+1,p));
                    end
                end
            end
        end
        
        % Check energy balance for pipes
        for p = 1:PipeCount
            node1 = LinkFromTo(p,1);
            node2 = LinkFromTo(p,2);
            
            if node1 <= JunctionCount && node2 <= JunctionCount
                % Calculate head difference
                dH = TrueState.JunctionHead(t+1,node1) - TrueState.JunctionHead(t+1,node2);
                energy_balance_errors(t,p) = abs(dH);
            end
        end
    end
    
    % Plot results
    figure('Position', [100 100 1200 400]);
    
    % Mass Balance Plot
    subplot(1,2,1);
    plot(1:24, mean(mass_balance_errors,2), 'k-', 'LineWidth', 2);
    xlabel('Time (hours)');
    ylabel('Average Mass Balance Error (GPM)');
    title('Mass Balance Validation');
    grid on;
    
    % Also print statistics
    fprintf('\nMass Balance Statistics:\n');
    fprintf('Maximum Error: %.6f GPM\n', max(mass_balance_errors(:)));
    fprintf('Average Error: %.6f GPM\n', mean(mass_balance_errors(:)));
    fprintf('Number of nodes with error > 1 GPM: %d\n', sum(max(mass_balance_errors) > 1));
    
    % Find worst junctions
    [max_errors, worst_junctions] = max(mass_balance_errors);
    [sorted_errors, time_indices] = sort(max_errors, 'descend');
    fprintf('\nTop 5 worst junction imbalances:\n');
    for i = 1:min(5, length(sorted_errors))
        t_idx = time_indices(i);
        j_idx = worst_junctions(t_idx);
        fprintf('Junction %d at time %d: %.6f GPM\n', j_idx, t_idx, sorted_errors(i));
    end
    
    % Energy Balance Plot
    subplot(1,2,2);
    plot(1:24, mean(energy_balance_errors,2), 'k-', 'LineWidth', 2);
    xlabel('Time (hours)');
    ylabel('Average Energy Balance Error (ft)');
    title('Energy Balance Validation');
    grid on;
    
    % Save detailed results
    detailed_results = struct();
    detailed_results.mass_balance_errors = mass_balance_errors;
    detailed_results.energy_balance_errors = energy_balance_errors;
    detailed_results.time_series = 1:24;
    save('network_balance_check.mat', 'detailed_results');
end