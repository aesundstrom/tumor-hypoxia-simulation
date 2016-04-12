function [] = simulator()

%% load operating parameters

load_operating_parameters;

%% declare parameters and data structures

% define essential simulation arrays
global occupation concentration impact_factor_cond replaceable_cond reproductive_cond fitness_1 fitness_n;
occupation          = zeros(x_dim, y_dim, z_dim);                      % cell types
concentration       = zeros(x_dim, y_dim, z_dim, num_particle_types);  % particle concentrations
impact_factor_cond  = zeros(x_dim, y_dim, z_dim, num_particle_types);  % impact factor of each particle type for each cell
replaceable_cond    = zeros(x_dim, y_dim, z_dim);                      % replacement status of each cell
reproductive_cond   = zeros(x_dim, y_dim, z_dim);                      % reproductive status of each cell
fitness_1           = zeros(x_dim, y_dim, z_dim);                      % local fitness (used in conjunction with occupation; each coordinate has one local fitness value)
fitness_n           = zeros(x_dim, y_dim, z_dim, num_cell_types);      % neighborhood fitness (each coordinate has a neighborhood fitness value for each cell type)

% define statistics arrays
global population fitness_1_avg fitness_1_std fitness_n_avg fitness_n_std diameter_x diameter_y diameter_z epc wallclock;
population          = zeros(num_cell_types, num_iters);  % population of each cell type
fitness_1_avg       = zeros(num_cell_types, num_iters);  % mean local fitness for each cell type
fitness_1_std       = zeros(num_cell_types, num_iters);  % standard deviation local fitness for each cell type
fitness_n_avg       = zeros(num_cell_types, num_iters);  % mean neighborhood fitness for each cell type
fitness_n_std       = zeros(num_cell_types, num_iters);  % standard deviation neighborhood fitness for each cell type
diameter_x          = zeros(num_cell_types, num_iters);  % global x-extent for each cell type
diameter_y          = zeros(num_cell_types, num_iters);  % global y-extent for each cell type
diameter_z          = zeros(num_cell_types, num_iters);  % global z-extent for each cell type
epc                 = zeros(num_cell_types, num_iters);  % global Euler-Poincare characteristic for each cell type
wallclock           = zeros(num_iters);                  % wallclock time for each tick's processing

%% set initial conditions

% define some useful landmarks
global mid_x mid_y mid_z mid_mid_x mid_mid_y mid_mid_z;
mid_x     = floor(x_dim / 2);
mid_y     = floor(y_dim / 2);
mid_z     = floor(z_dim / 2) + 1;
mid_mid_x = floor(x_dim / 4);
mid_mid_y = floor(y_dim / 4);
mid_mid_z = floor(z_dim / 4) + 1;

% initialize particle concentrations by particle type
for pt = 1 : num_particle_types
    concentration(:,:,:,pt) = init_concentration(pt);
end

% initialize the cell type occupations
occupation_delay_elapsed = false;
if ~delay_occupation_by
    initialize_occupation;
    occupation_delay_elapsed = true;
end

%% setup for outputting results

if output_results

    % create timestamp
    [year, month, day, hour, minute, second] = datevec(now);
    years   = num2str(year);
    months  = num2str(month);
    days    = num2str(day);
    hours   = num2str(hour);
    minutes = num2str(minute);
    seconds = num2str(floor(second));
    if length(months) == 1
        months = strcat('0', months);
    end
    if length(days) == 1
        days = strcat('0', days);
    end
    if length(hours) == 1
        hours = strcat('0', hours);
    end
    if length(minutes) == 1
        minutes = strcat('0', minutes);
    end
    if length(seconds) == 1
        seconds = strcat('0', seconds);
    end
    ts = sprintf('%s_%s_%s__%s_%s_%s', years, months, days, hours, minutes, seconds);
    
    % create new results directory and copy this code into it
    global path;
    path = strcat('/tmp/simulations/', ts);
    mkdir(path);
    copyfile('simulate_metabolic_symbiosis_2d.m', path);

end

%% simulation loop

for tick = 1 : num_iters
    
    %% prologue
    
    % set the timer for this tick
    tic;
    
    % set the time scale separation parameter
    tau = tick * d_tau;
    
    % initialize once the cell type occupations after the delay
    if tick > delay_occupation_by && ~occupation_delay_elapsed
        initialize_occupation;
        occupation_delay_elapsed = true;
    end
    
    % initialize the conditional matrices after the delay
    if tick > delay_occupation_by
        impact_factor_cond = reshape(impact_factor(occupation,:), x_dim, y_dim, z_dim, num_particle_types);
        replaceable_cond   = replaceable(occupation);
        reproductive_cond  = reproductive(occupation);
    end

    %% phase I: consume and release particles
    
    % if past the occupation delay
    if tick > delay_occupation_by
        
        % DEFAULT: consume and release particles according to the cells' respective consumption and release rates
        for pt = 1 : num_particle_types
            concentration_pt = concentration(:,:,:,pt);
            for ct = 1 : num_cell_types
                concentration_pt(occupation == ct) = concentration_pt(occupation == ct) * (1 - consume_rate(ct,pt) + release_rate(ct,pt));
            end
            concentration(:,:,:,pt) = concentration_pt;
        end
        
        % CONDITIONAL: for each cell type defined in the conditionals
        for cond_key = conditionals.keys
            
            % define the cell type from the conditional key
            ct = cond_key{1};
            
            % fetch the conditional block
            cond = conditionals(ct);
            
            % for each trigger-action defined for this cell type
            for ta = 1 : numel(cond)
                
                % parse out the trigger set and the action set
                triggers = cond{ta}{1};  % trigger set
                actions  = cond{ta}{2};  % action set
                
                % declare a result matrix for this trigger set
                ts_result = ones(x_dim, y_dim, z_dim);
                
                % for each trigger in the trigger set
                for t = 1 : numel(triggers)
                    
                    % fetch the trigger out of the trigger set
                    trigger = triggers{t};
                    
                    % parse out the trigger elements
                    trigger_pt = trigger{1};  % this trigger's particle type
                    trigger_op = trigger{2};  % this trigger's operator
                    trigger_pc = trigger{3};  % this trigger's particle concentration
                    
                    % obtain the appropriate concentrations for this trigger's particle type
                    concentration_trigger_pt = concentration(:,:,:,trigger_pt);
                    
                    % determine the results of this trigger based on this trigger's operator
                    switch trigger_op
                        case '<'
                            ts_result = ts_result & (occupation == ct) & (concentration_trigger_pt <  trigger_pc);
                        case '>'
                            ts_result = ts_result & (occupation == ct) & (concentration_trigger_pt >  trigger_pc);
                        case '='
                            ts_result = ts_result & (occupation == ct) & (concentration_trigger_pt == trigger_pc);
                    end
                    
                end  % for each trigger in the trigger set
                
                % for each action defined for this trigger set
                for a = 1 : numel(actions)
                    
                    % fetch the action out of the action set
                    action = actions{a};
                    
                    % parse out the action elements
                    action_op = action{1};  % this action's operator
                    if ~isequal(action_op, 'a')
                        action_o1 = action{2};  % this action's operand 1 (e.g., boolean target value, particle type, jump cell type)
                    end
                    if ~isequal(action_op, 'a') && ~isequal(action_op, 'j') && ~isequal(action_op, '-') && ~isequal(action_op, '+')
                        action_o2 = action{3};  % this action's operand 2 (e.g., particle concentration)
                    end
                    
                    % for the trigger set result cells, perform the action
                    switch action_op
                        case 'a'
                            occupation(ts_result) = 2;
                            for pt = 1 : num_particle_types
                                impact_factor_cond_pt = impact_factor_cond(:,:,:,pt);
                                impact_factor_cond_pt(ts_result) = impact_factor(2,pt);
                                impact_factor_cond(:,:,:,pt) = impact_factor_cond_pt;
                            end
                            replaceable_cond(ts_result) = replaceable(2);
                            reproductive_cond(ts_result) = reproductive(2);
                        case 'j'
                            occupation(ts_result) = action_o1;
                            for pt = 1 : num_particle_types
                                impact_factor_cond_pt = impact_factor_cond(:,:,:,pt);
                                impact_factor_cond_pt(ts_result) = impact_factor(action_o1,pt);
                                impact_factor_cond(:,:,:,pt) = impact_factor_cond_pt;
                            end
                            replaceable_cond(ts_result) = replaceable(action_o1);
                            reproductive_cond(ts_result) = reproductive(action_o1);
                        case '-'
                            replaceable_cond(ts_result) = action_o1;
                        case '+'
                            reproductive_cond(ts_result) = action_o1;
                        case 'c'
                            concentration_o1 = concentration(:,:,:,action_o1);
                            concentration_o1(ts_result) = concentration_o1(ts_result) * (1 - action_o2);
                            concentration(:,:,:,action_o1) = concentration_o1;
                        case 'r'
                            concentration_o1 = concentration(:,:,:,action_o1);
                            concentration_o1(ts_result) = concentration_o1(ts_result) * (1 + action_o2);
                            concentration(:,:,:,action_o1) = concentration_o1;
                        case 'i'
                            impact_factor_cond_o1 = impact_factor_cond(:,:,:,action_o1);
                            impact_factor_cond_o1(ts_result) = action_o2;
                            impact_factor_cond(:,:,:,action_o1) = impact_factor_cond_o1;
                    end
                    
                end  % for each action defined for this trigger set
                
            end  % for each trigger-action defined for this cell type
            
        end  % for each cell type defined in the conditionals
        
    end  % if tick > delay_occupation
        
    %% phase II: diffuse particles
    
    % diffuse particles according to the particles' respective diffusion rates
    for pt = 1 : num_particle_types
        sigma = sqrt(2 * diffusion_rate(pt));
        if z_dim > 1
            G = gauss3([5 5 5], [sigma sigma sigma]);
        else
            G = fspecial('gaussian', [5 5], sigma);
        end
        concentration(:,:,:,pt) = imfilter(concentration(:,:,:,pt), G, 'replicate', 'same', 'conv');
    end
    
    % restore concentrations to within their basal bounds
    for pt = 1 : num_particle_types
        concentration_pt = concentration(:,:,:,pt);
        for ct = 1 : num_cell_types
            bl = basal_lower(ct, pt);
            bu = basal_upper(ct, pt);
            concentration_pt((occupation == ct) & (concentration_pt < bl)) = bl;
            concentration_pt((occupation == ct) & (concentration_pt > bu)) = bu;
        end
        concentration(:,:,:,pt) = concentration_pt;
    end

    %% phase III: compute cell fitness
    
    % if past the occupation delay
    if tick > delay_occupation_by
        
        % compute local fitness of each voxel
        fitness_1 = zeros(x_dim, y_dim, z_dim);
        for i = 1 : x_dim
            for j = 1 : y_dim
                for k = 1 : z_dim
                    
                    % get the cell type at this voxel
                    cell_type = occupation(i,j,k);
                    
                    % skip vessel and empty cell types, since these cannot obtain a local fitness
                    if cell_type == 1 || cell_type == 2
                        continue;
                    end
                    
                    % comute local fitness based on particle type concentrations and their respective impact factors
                    impact_factors = reshape(impact_factor_cond(i,j,k,:), 1, num_particle_types);
                    concentrations = reshape(concentration(i,j,k,:), num_particle_types, 1);
                    fitness_1(i,j,k) = impact_factors * concentrations;
                    
                end
            end
        end
        
        % compute neighborhood fitness of each voxel
        fitness_n = zeros(x_dim, y_dim, z_dim, num_cell_types);
        for i = 1 : x_dim
            for j = 1 : y_dim
                for k = 1 : z_dim
                    
                    % get the cell type at this voxel
                    cell_type = occupation(i,j,k);
                    
                    % skip vessel cell type, since they cannot be targeted for a new cell type (but empty cells still can)
                    if cell_type == 1
                        continue;
                    end
                    
                    % compute neighborhood fitness; loop from cell type 3 (post-vessel, post-empty) to the last type
                    for ct = 3 : num_cell_types
                        num_neighbors = 0;
                        sum_neighbors = 0;
                        for ii = i-1 : i+1
                            if ii < 1 || ii > x_dim
                                continue;
                            end
                            for jj = j-1 : j+1
                                if jj < 1 || jj > y_dim
                                    continue;
                                end
                                for kk = k-1 : k+1
                                    if kk < 1 || kk > z_dim
                                        continue;
                                    end
                                    if occupation(ii,jj,kk) ~= ct
                                        continue;
                                    end
                                    num_neighbors = num_neighbors + 1;
                                    sum_neighbors = sum_neighbors + fitness_1(ii,jj,kk);
                                end  % for kk
                            end  % for jj
                        end  % for ii
                        if num_neighbors == 0
                            fitness_n(i,j,k,ct) = 0;
                        else
                            fitness_n(i,j,k,ct) = sum_neighbors / num_neighbors;
                        end
                    end  % for ct
                    
                end  % for k
            end  % for j
        end  % for i
        
    end  % if tick > delay_occupation
    
    %% interlude: compute and plot statistics
    
    % compute statistics
    for ct = 1 : num_cell_types
        
        % occupations for this cell type
        occupation_ct = occupation == ct;
        
        % population of this cell type
        population(ct, tick) = sum(sum(sum(occupation_ct)));
        
        % local fitness of this cell type
        fitness_1_avg(ct, tick) = mean(fitness_1(occupation_ct));
        fitness_1_std(ct, tick) = std(fitness_1(occupation_ct));
        
        % neighborhood fitness of this cell type
        fitness_n_ct            = fitness_n(:,:,:,ct);
        fitness_n_avg(ct, tick) = mean(fitness_n_ct(occupation > 1));
        fitness_n_std(ct, tick) = std(fitness_n_ct(occupation > 1));
        
        % x,y,z-extent for this cell type
        [cti, ctj, ctk] = ind2sub(size(occupation_ct), find(occupation_ct)); 
        if numel(cti) == 0
            diameter_x(ct, tick) = 0;
        else
            diameter_x(ct, tick) = max(cti) - min(cti) + 1;
        end
        if numel(ctj) == 0
            diameter_y(ct, tick) = 0;
        else
            diameter_y(ct, tick) = max(ctj) - min(ctj) + 1;
        end
        if numel(ctk) == 0
            diameter_z(ct, tick) = 0;
        else
            diameter_z(ct, tick) = max(ctk) - min(ctk) + 1;
        end
        
        % Euler-Poincare characteristic for this cell type
        if z_dim > 1
            epc(ct, tick) = imEuler3d(occupation_ct);
        else
            epc(ct, tick) = imEuler2d(occupation_ct);
        end
        
    end
  
    % plot statistics
    if ~mod(tick, plot_every)
        plot_statistics();
    end

    %% phase IV: reproduce cells probabilistically
    
    % if past the occupation delay
    if tick > delay_occupation_by
        
        % if it's time to reproduce
        if ~mod(tick, reproduce_every)
            
            % probabilistically determine the fate of each voxel
            for i = 1 : x_dim
                for j = 1 : y_dim
                    for k = 1 : z_dim
                        
                        % if this cell type is not replaceable, then skip it
                        if ~replaceable_cond(i,j,k)
                            continue;
                        end
                        
                        % declare a probability vector; the first two elements (representing vessel and empty cell types) won't be used
                        probabilities = zeros(1, num_cell_types);
                        
                        % obtain the probabilities from the neighborhood fitness values; loop from cell type 3 (post-vessel, post-empty) to the last type
                        for ct = 3 : num_cell_types
                            probabilities(ct) = fitness_n(i,j,k,ct);
                        end
                        
                        % sum the probabilities
                        sum_probabilities = sum(probabilities);

                        % if the sum of probabilities is greater than one, then normalize the probabilities to one
                        shrinkage_factor = 1;
                        if sum_probabilities > 1
                            shrinkage_factor = 1 / sum_probabilities;
                        end
                        probabilities = shrinkage_factor * probabilities;
                        
                        % partition the [0,1] interval by the set of probabilities, contiguously placed along the interval, and randomly point into the interval to select the target cell type for this cell; loop from cell type 3 (post-vessel, post-empty) to the last type
                        essay = rand;
                        beg_num = 0;
                        end_num = 0;
                        target_cell_type = 0;
                        for ct = 3 : num_cell_types
                            beg_num = end_num;
                            end_num = end_num + probabilities(ct);
                            if essay >= beg_num && essay < end_num
                                target_cell_type = ct;
                                break;
                            end
                        end
                        
                        % if the random point lies beyond the last cell type's probability range, then it must become an empty cell
                        if target_cell_type == 0
                            target_cell_type = 2;
                        end
                        
                        % if any of the neighboring cells are of the target type and are reproductive, then mutate the current cell type to the target cell type
                        neighbor_is_of_target_cell_type_and_reproductive = false;
                        for ii = i-1 : i+1
                            if ii < 1 || ii > x_dim
                                continue;
                            end
                            for jj = j-1 : j+1
                                if jj < 1 || jj > y_dim
                                    continue;
                                end
                                for kk = k-1 : k+1
                                    if kk < 1 || kk > z_dim
                                        continue;
                                    end
                                    if (occupation(ii,jj,kk) == target_cell_type) && reproductive_cond(ii,jj,kk)
                                        neighbor_is_of_target_cell_type_and_reproductive = true;
                                    end
                                end  % for kk
                            end  % for jj
                        end  % for ii
                        if neighbor_is_of_target_cell_type_and_reproductive
                            occupation(i,j,k) = target_cell_type;
                        end
                        
                    end  % for k
                end  % for j
            end  % for i
            
        end  % reproduce every
    
    end  % if tick > delay_occupation
    
    %% epilogue
    
    % record the timer for this tick
    wallclock(tick) = toc;
    
end  % for tick

end  % function simulator
