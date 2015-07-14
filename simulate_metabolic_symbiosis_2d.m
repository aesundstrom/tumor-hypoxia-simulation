function [] = simulate_metabolic_symbiosis_2d()

%% declare parameters and data structures

% global parameters
s_dim               = 40;     % symmetric dimension (used for convenience, when x, y, z dims are identical)
x_dim               = s_dim;  % x dimension
y_dim               = s_dim;  % y dimension
z_dim               = 1;      % z dimension
num_cell_types      = 4;      % number of cell types
num_particle_types  = 3;      % number of particle types
num_iters           = 1000;   % number of simulation clock ticks
d_tau               = 1;      % time scale separation paremeter
plot_every          = 1;      % number of clock ticks between plottings
reproduce_every     = 1;      % number of clock ticks between probablilistic reproductions
delay_occupation_by = 0;      % number of clock ticks before initializing cell type occupations (used to establish gradients prior to exposing cells)
plot_3d             = 1;      % predicate for plotting 3D occupation per cell type

% output predicates
output_results                = 1;        % predicate for outputting plots (global switch for those below)
output_cell_types             = 1;        % predicate for outputting cell types
output_time_series            = 1;        % predicate for outputting time series
output_local_fitness_3d       = 1;        % predicate for outputting local fitness in 3D
output_local_fitness          = 1;        % predicate for outputting local fitness
output_neighborhood_fitness   = 0;        % predicate for outputting neighborhood fitness
output_particle_concentration = 0;        % predicate for outputting particle concentration
output_only_cell_types        = [3 4];    % subset of cell types to output
output_only_particle_types    = [1 2 3];  % subset of particle types to output

% diffusion rate of each particle type
diffusion_rate      = [1    1    1];

% initial concentration of each particle type
init_concentration  = [0    0    0.01];

% basal lower bound of each particle type for each cell type
basal_lower         = [0.1  0.1  0;   % vessel
                       0.1  0.1  0;   % empty
                       0.1  0.1  0;   % alpha
                       0.1  0.1  0];  % beta

% basal upper bound of each particle type for each cell type
basal_upper         = [Inf  Inf  Inf;   % vessel
                       Inf  Inf  Inf;   % empty
                       Inf  Inf  Inf;   % alpha
                       Inf  Inf  Inf];  % beta
                       
% consume rate of each particle type for each cell type
consume_rate        = [0    0    10.0;  % vessel
                       0    0    0;     % empty
                       1.0  0    0;     % alpha
                       0    1.0  1.0];  % beta

% release rate of each particle type for each cell type
release_rate        = [10.0 10.0 0;    % vessel
                       0    0    0;    % empty
                       0    0    1.0;  % alpha
                       0    0    0];   % beta
                       
% impact factor of each particle type for each cell type
impact_factor       = [0    0    0;   % vessel
                       0    0    0;   % empty
                       1    0    0;   % alpha
                       0    1    1];  % beta

% replacement status of each cell type
replaceable         = [0;   % vessel
                       1;   % empty
                       1;   % alpha
                       1];  % beta

% reproductive status of each cell type
reproductive        = [0;   % vessel
                       1;   % empty
                       1;   % alpha
                       1];  % beta
                   
% color of each cell type
cell_type_color_map = [1    1    1;   % vessel (white)
                       0    0    1;   % empty  (blue)
                       1    0    0;   % alpha  (red)
                       0    1    0];  % beta   (green)

% declare conditional triggers and actions:
%     trigger = { <particle type> <comparative operator = {<,=,>}> <value> }
%     action  = { <operator = {a,j,-,+,c,r,i}> <operand 1 = {jump cell type, boolean target value, particle type}> <operand 2 = {particle concentration, impact factor}> }
% note:
%     operator:  number of operands:
%     {a}        0
%     {j,-,+}    1
%     {c,r,i}    2
% example conditional trigers and actions:
%     conditionals(n) = { { { {t1} {t2} {t3} } { {a1} {a2}      } }
%                         { { {t1} {t2}      } { {a1} {a2} {a3} } }
%                       };
conditionals    = containers.Map('KeyType', 'int32', 'ValueType', 'any');
              
% define essential simulation arrays 
occupation          = zeros(x_dim, y_dim, z_dim);                      % cell types
concentration       = zeros(x_dim, y_dim, z_dim, num_particle_types);  % particle concentrations
impact_factor_cond  = zeros(x_dim, y_dim, z_dim, num_particle_types);  % impact factor of each particle type for each cell
replaceable_cond    = zeros(x_dim, y_dim, z_dim);                      % replacement status of each cell
reproductive_cond   = zeros(x_dim, y_dim, z_dim);                      % reproductive status of each cell
fitness_1           = zeros(x_dim, y_dim, z_dim);                      % local fitness (used in conjunction with occupation; each coordinate has one local fitness value)
fitness_n           = zeros(x_dim, y_dim, z_dim, num_cell_types);      % neighborhood fitness (each coordinate has a neighborhood fitness value for each cell type)

% define statistics arrays
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
    path = strcat('/tmp/simulations/', ts);
    mkdir(path);
    copyfile('gd.m', path);
    
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

    %% plot the statistics of interest
    function [] = plot_statistics()
        
        % define some useful constants
        gray_256    = grade([1 1 1], 256);
        num_rows    = 6;
        num_ct_cols = (num_cell_types - 2) + 1;  % no columns for vessel and empty cell types + 1 column for corresponding time series
        num_cols    = max(num_ct_cols, num_particle_types);
        del_cols    = num_cols - num_ct_cols;
        ct_range    = 3 : num_cell_types;
        pt_range    = 1 : num_particle_types;
        cursor      = 1;
        
        % set up the figure
        figure(1);
        clf;
        
        % cell types (2-D slice)
        subplot(num_rows, num_cols, cursor);
        imagesc(occupation(:,:,mid_z), [1 size(cell_type_color_map,1)]);
        colormap(cell_type_color_map);
        freezeColors;
        set(gca, 'YDir', 'normal');
        axis square;
        cursor = cursor + 1;
        
        % performance (time series)
        subplot(num_rows, num_cols, cursor);
        domain = 1:tick-1;
        wc     = wallclock(domain);
        wc_avg = mean(wc);
        wc_std = std(wc);
        plot(domain, wc);
        freezeColors;
        title(['avg=' sprintf('%3.3f', wc_avg) ', std=' sprintf('%3.3f', wc_std)]);
        axis square;
        cursor = cursor + 1;
        
        % pad, if necessary
        pad_cursor();
        
        % cell type 3D plot
        for ct = ct_range
            if plot_3d
                subplot(num_rows, num_cols, cursor);
                occ_ct = occupation == ct;
                fit_ct = fitness_1(occ_ct);
                min_fit_ct = min(fit_ct);
                max_fit_ct = max(fit_ct);
                if min_fit_ct < max_fit_ct
                    norm_fit_ct = (fit_ct - min_fit_ct) ./ (max_fit_ct - min_fit_ct);
                else
                    norm_fit_ct = 0;
                end
                color_density = 256;
                color_grade = grade(cell_type_color_map(ct,:), color_density);
                color_idx = floor((color_density - 1) * norm_fit_ct) + 1;
                scatter_colors = color_grade(color_idx,:);
                [x,y,z] = ind2sub(size(occ_ct), find(occ_ct));
                % note the axis order for plotting: <y,x,z>
                scatter3(y, x, z, 20, scatter_colors, 'filled', 's');
                colormap(grade(cell_type_color_map(ct,:), color_density));
                ch = colorbar;
                if min_fit_ct < max_fit_ct
                    caxis([min_fit_ct max_fit_ct]);
                end
                cbfreeze(ch);
                xlabel('x');
                ylabel('y');
                zlabel('z');
                if z_dim > 1
                    axis([1 x_dim 1 y_dim 1 z_dim]);
                else
                    axis([1 x_dim 1 y_dim]);
                end
                axis square;
            end
            cursor = cursor + 1;
        end
        
        % population by cell type (time series)
        subplot(num_rows, num_cols, cursor);
        hold on;
        for ct = 1 : num_cell_types
            plot(1:tick, population(ct, 1:tick), 'color', cell_type_color_map(ct,:));
        end
        hold off;
        freezeColors;
        if tick > 1
            xlim([1 tick]);
        end
        axis square;
        cursor = cursor + 1;
                
        % pad, if necessary
        pad_cursor();
        
        % cell type local fitness (2-D slice)
        for ct = ct_range
            subplot(num_rows, num_cols, cursor);
            imagesc((occupation(:,:,mid_z) == ct) .* fitness_1(:,:,mid_z));
            colormap(grade(cell_type_color_map(ct,:), 256));
            cbfreeze(colorbar);
            freezeColors;
            set(gca, 'YDir', 'normal');
            axis square;
            cursor = cursor + 1;
        end
        
        % local fitness by cell type (time series)
        subplot(num_rows, num_cols, cursor);
        plot_time_series(fitness_1_avg(ct_range,:), fitness_1_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
        axis square;
        cursor = cursor + 1;
        
        % pad, if necessary
        pad_cursor();

        % cell type neighborhood fitness (2-D slice)
        for ct = ct_range
            subplot(num_rows, num_cols, cursor);
            imagesc((occupation(:,:,mid_z) > 1) .* fitness_n(:,:,mid_z,ct));
            colormap(grade(cell_type_color_map(ct,:), 256));
            cbfreeze(colorbar);
            freezeColors;
            set(gca, 'YDir', 'normal');
            axis square;
            cursor = cursor + 1;
        end
                
        % neighborhood fitness by cell type (time series)
        subplot(num_rows, num_cols, cursor);        
        plot_time_series(fitness_n_avg(ct_range,:), fitness_n_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
        axis square;
        cursor = cursor + 1;
        
        % pad, if necessary
        pad_cursor();
        
        % cell type x,y,z-extent
        for ct = ct_range
            subplot(num_rows, num_cols, cursor);
            hold on;
            % note the axis order for plotting: <y,x,z>
            plot(1:tick, diameter_y(ct, 1:tick), 'r');
            plot(1:tick, diameter_x(ct, 1:tick), 'g');
            plot(1:tick, diameter_z(ct, 1:tick), 'b');
            hold off;
            if tick > 1
                xlim([1 tick]);
            end
            axis square;
            cursor = cursor + 1;
        end
        
        % Euler-Poincare characteristic by cell type (time series)
        subplot(num_rows, num_cols, cursor);
        hold on;
        for ct = 1 : num_cell_types
            plot(1:tick, epc(ct, 1:tick), 'color', cell_type_color_map(ct,:));
        end
        hold off;
        freezeColors;
        axis square;
        cursor = cursor + 1;
        
        % pad, if necessary
        pad_cursor();

        % particle type concentrations (2-D slice)
        for pt = pt_range
            subplot(num_rows, num_cols, cursor);
            con_pt = concentration(:,:,mid_z,pt);
            min_con_pt = min(min(con_pt));
            max_con_pt = max(max(con_pt));
            mesh(con_pt);
            colormap(gray_256);
            min_rel = '=';
            if min_con_pt > 1000
                min_rel = '>';
                min_con_pt = 1000;
            end
            max_rel = '=';
            if max_con_pt > 1000
                max_rel = '>';
                max_con_pt = 1000;
            end
            title(['min' min_rel sprintf('%3.3f', min_con_pt) ', max' max_rel sprintf('%3.3f', max_con_pt)]);
            axis square;
            xlim([1 x_dim]);
            ylim([1 y_dim]);
            cursor = cursor + 1;
        end
        
        % if outputting results, then create and save a stand-alone figure
        if output_results
            if output_cell_types
                save_cell_types();
            end
            if output_time_series
                save_time_series(ct_range);
            end
            if output_local_fitness_3d
                save_local_fitness_3d(ct_range);
            end
            if output_local_fitness
                save_local_fitness(ct_range);
            end
            if output_neighborhood_fitness
                save_neighborhood_fitness(ct_range);
            end
            if output_particle_concentration
                save_particle_concentration(gray_256, pt_range)
            end
        end
        
        function [rc] = row_cursor(cursor, num_cols)
            rc = mod(cursor, num_cols);
            if rc == 0
                rc = num_cols;
            end
        end  % function row_cursor
        
        function [] = pad_cursor()
            for p = 1 : (num_cols - row_cursor(cursor - 1, num_cols))
                cursor = cursor + 1;
            end
        end  % function pad_cursor
        
    end  % function plot_statistics

    %% plot cell types and save them to disk
    function [] = save_cell_types()
        
        % set up an invisible figure
        h = figure('Visible', 'off');

        % cell types (2-D slice)
        imagesc(occupation(:,:,mid_z), [1 size(cell_type_color_map,1)]);
        colormap(cell_type_color_map);
        freezeColors;
        set(gca, 'YDir', 'normal');
        axis square;
        title({'Cell Population' ['t=' num2str(tick)]});
        
        % create figure file name and save figure
        len_num_iters = length(num2str(num_iters));
        len_tick = length(num2str(tick));
        len_pad = len_num_iters - len_tick;
        pad = '';
        for pi = 1 : len_pad
            pad = strcat('0', pad);
        end
        fig_file = strcat(path, '/', 'cells_', pad, num2str(tick));
        saveas(gcf, fig_file, 'pdf');
        
        % delete figure
        close(h);
        
    end  % function save_cell_types

    %% plot certain time series and save them to disk
    function [] = save_time_series(ct_range)
        
        % set up an invisible figure
        h = figure('Visible', 'off');

        % population by cell type (time series)
        subplot(3,1,1);
        hold on;
        for ct = 1 : num_cell_types
            plot(1:tick, population(ct, 1:tick), 'color', cell_type_color_map(ct,:));
        end
        hold off;
        if tick > 1
            xlim([1 tick]);
        end
        title('Time Evolution of Populations by Cell Type');
        
        % local fitness by cell type (time series)
        subplot(3,1,2);
        plot_time_series(fitness_1_avg(ct_range,:), fitness_1_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
        title('Time Evolution of Local Fitness by Cell Type');
        
        % neighborhood fitness by cell type (time series)
        subplot(3,1,3);
        plot_time_series(fitness_n_avg(ct_range,:), fitness_n_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
        title('Time Evolution of Neighborhood Fitness by Cell Type');
        
        % create figure file name and save figure
        fig_file = strcat(path, '/', 'time_series');
        saveas(h, fig_file, 'fig');
        saveas(h, fig_file, 'pdf');
        
        % delete figure
        close(h);
        
    end  % function save_time_series

    %% plot local fitness 3D and save them to disk
    function [] = save_local_fitness_3d(ct_range)
        
        % cell type 3D plot
        for ct = ct_range
            
            % process only designated cell types
            if ~numel(find(output_only_cell_types == ct))
                continue;
            end
            
            % set up an invisible figure
            h = figure('Visible', 'off');
            
            % plot the image
            occ_ct = occupation == ct;
            fit_ct = fitness_1(occ_ct);
            min_fit_ct = min(fit_ct);
            max_fit_ct = max(fit_ct);
            norm_fit_ct = fit_ct;
            if min_fit_ct < max_fit_ct
                norm_fit_ct = (fit_ct - min_fit_ct) ./ (max_fit_ct - min_fit_ct);
            end
            color_density = 256;
            color_grade = grade(cell_type_color_map(ct,:), color_density);
            color_idx = floor((color_density - 1) * norm_fit_ct) + 1;
            scatter_colors = color_grade(color_idx,:);
            [x,y,z] = ind2sub(size(occ_ct), find(occ_ct));
            % note the axis order for plotting: <y,x,z>
            scatter3(y, x, z, 20, scatter_colors, 'filled', 's');
            colormap(grade(cell_type_color_map(ct,:), color_density));
            drawnow;
            ch = colorbar;
            if min_fit_ct < max_fit_ct
                caxis([min_fit_ct max_fit_ct]);
            end
            cbfreeze(ch);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            if z_dim > 1
                axis([1 x_dim 1 y_dim 1 z_dim]);
            else
                axis([1 x_dim 1 y_dim]);
            end
            axis square;
            title({['Local Fitness of Cell Type ' num2str(ct) ' (3D View)'] ['t=' num2str(tick)]});
            
            % create figure file name and save figure
            len_num_iters = length(num2str(num_iters));
            len_tick = length(num2str(tick));
            len_pad = len_num_iters - len_tick;
            pad = '';
            for pi = 1 : len_pad
                pad = strcat('0', pad);
            end
            fig_file = strcat(path, '/', 'local_fitness_3d_', num2str(ct), '_', pad, num2str(tick));
            saveas(gcf, fig_file, 'pdf');
            
            % delete figure
            close(h);
            
        end
                
    end  % function save_local_fitness_3d

    %% plot local fitness and save them to disk
    function [] = save_local_fitness(ct_range)
        
        % cell type local fitness (2-D slice)
        for ct = ct_range
            
            % process only designated cell types
            if ~numel(find(output_only_cell_types == ct))
                continue;
            end
            
            % set up an invisible figure
            h = figure('Visible', 'off');

            % plot the image
            imagesc((occupation(:,:,mid_z) == ct) .* fitness_1(:,:,mid_z));
            colormap(grade(cell_type_color_map(ct,:), 256));
            cbfreeze(colorbar);
            freezeColors;
            set(gca, 'YDir', 'normal');
            axis square;
            title({['Local Fitness of Cell Type ' num2str(ct)] ['t=' num2str(tick)]});

            % create figure file name and save figure
            len_num_iters = length(num2str(num_iters));
            len_tick = length(num2str(tick));
            len_pad = len_num_iters - len_tick;
            pad = '';
            for pi = 1 : len_pad
                pad = strcat('0', pad);
            end
            fig_file = strcat(path, '/', 'local_fitness_', num2str(ct), '_', pad, num2str(tick));
            saveas(gcf, fig_file, 'pdf');
        
            % delete figure
            close(h);

        end
        
    end  % function save_local_fitness

    %% plot neighborhood fitness and save them to disk
    function [] = save_neighborhood_fitness(ct_range)
        
        % cell type neighborhood fitness (2-D slice)
        for ct = ct_range
            
            % process only designated cell types
            if ~numel(find(output_only_cell_types == ct))
                continue;
            end
            
            % set up an invisible figure
            h = figure('Visible', 'off');

            % plot the image
            imagesc((occupation(:,:,mid_z) > 1) .* fitness_n(:,:,mid_z,ct));
            colormap(grade(cell_type_color_map(ct,:), 256));
            cbfreeze(colorbar);
            freezeColors;
            set(gca, 'YDir', 'normal');
            axis square;
            title({['Neighborhood Fitness of Cell Type ' num2str(ct)] ['t=' num2str(tick)]});

            % create figure file name and save figure
            len_num_iters = length(num2str(num_iters));
            len_tick = length(num2str(tick));
            len_pad = len_num_iters - len_tick;
            pad = '';
            for pi = 1 : len_pad
                pad = strcat('0', pad);
            end
            fig_file = strcat(path, '/', 'neighborhood_fitness_', num2str(ct), '_', pad, num2str(tick));
            saveas(gcf, fig_file, 'pdf');
        
            % delete figure
            close(h);

        end
        
    end  % function save_neighborhood_fitness

    %% plot particle concentration and save them to disk
    function [] = save_particle_concentration(gray_256, pt_range)
        
        % particle type concentrations (2-D slice)
        for pt = pt_range
            
            % process only designated particle types
            if ~numel(find(output_only_particle_types == ct))
                continue;
            end
            
            % set up an invisible figure
            h = figure('Visible', 'off');

            % plot the image
            con_pt = concentration(:,:,mid_z,pt);
            min_con_pt = min(min(con_pt));
            max_con_pt = max(max(con_pt));
            mesh(con_pt);
            colormap(gray_256);
            min_rel = '=';
            if min_con_pt > 1000
                min_rel = '>';
                min_con_pt = 1000;
            end
            max_rel = '=';
            if max_con_pt > 1000
                max_rel = '>';
                max_con_pt = 1000;
            end
            axis square;
            xlim([1 x_dim]);
            ylim([1 y_dim]);
            title({['Concentration of Particle Type ' num2str(pt)] ['t=' num2str(tick)] ['min' min_rel sprintf('%3.3f', min_con_pt) ', max' max_rel sprintf('%3.3f', max_con_pt)]});
            
            % create figure file name and save figure
            len_num_iters = length(num2str(num_iters));
            len_tick = length(num2str(tick));
            len_pad = len_num_iters - len_tick;
            pad = '';
            for pi = 1 : len_pad
                pad = strcat('0', pad);
            end
            fig_file = strcat(path, '/', 'particle_concentration_', num2str(pt), '_', pad, num2str(tick));
            saveas(gcf, fig_file, 'pdf');
        
            % delete figure
            close(h);

        end
        
    end  % function save_particle_concentration

    %% plot a time series average +/- standard deviation (average curve surrounded by +/- gray patches)
    %  a: average            (m time series x n ticks)
    %  s: standard deviation (m time series x n ticks)
    %  c: color map          (m time series x 3 {r,g,b})
    %  b: begin time         (integer)
    %  e: end time           (integer)
    function [] = plot_time_series(a, s, c, b, e)
        domain = b : e;
        gray = [0.9 0.9 0.9];
        hold on;
        for t = 1 : size(a,1)
            patch([domain fliplr(domain)], [a(t,domain) - s(t,domain), fliplr(a(t,domain) + s(t,domain))], gray, 'LineStyle', 'none');
            plot(domain, a(t,domain), 'color', c(t,:));
        end
        hold off;
        freezeColors;
        if tick > 1
            xlim([1 tick]);
        end;
    end  % function plot_time_series

    %% create a color map using a basis and a granularity
    %  c: color map basis     ([r g b])
    %  n: density of gradient (integer)
    function [q] = grade(c, n)
        r = linspace(0,c(1),n)';
        g = linspace(0,c(2),n)';
        b = linspace(0,c(3),n)';
        q = horzcat(r,g,b);
    end  % function grade

    %% create a 3D Gaussian filter using a kernel and standard deviation
    %  k: kernel             ([kernel_dim_x kernel_dim_y kernel_dim_z])
    %  s: standard deviation ([sigma_x sigma_y sigma_z])
    function [g] = gauss3(k, s)
        k = floor(k/2);
        [gx, gy, gz] = ndgrid(-k(1):k(1), -k(2):k(2), -k(3):k(3));
        gx = gx / s(1);
        gy = gy / s(2);
        gz = gz / s(3);
        g = exp(-(gx .* gx + gy .* gy + gz .* gz) .* 0.5);
        g = g / sum(g(:));
    end  % function gauss3

    %% create a random 3D binary matrix using a gradient function for density along any combination of axis directions
    %  xd   : x dimension
    %  yd   : y dimension
    %  zd   : z dimension
    %  f    : gradient function (f(upper)...f(lower) guaranteed to be in the range [0,1])
    %  d_beg: domain lower bound
    %  d_end: domain upper bound
    %  dir  : direction of gradient (3-element binary vector)
    function [dm] = density_grad(xd, yd, zd, f, d_beg, d_end, dir)
        dm = [];
        
        domain = linspace(d_beg,d_end,xd);
        target = f(domain);
        dx = [];
        for z_idx = 1 : zd
            zm = [];
            for x_idx = 1 : xd
                zm = horzcat(zm, rand(yd,1) < target(x_idx));
            end
            dx(:,:,z_idx) = zm;
        end
                
        domain = linspace(d_beg,d_end,yd);
        target = f(domain);
        dy = [];
        for z_idx = 1 : zd
            zm = [];
            for y_idx = 1 : yd
                zm = vertcat(zm, rand(1,xd) < target(y_idx));
            end
            dy(:,:,z_idx) = zm;
        end
        
        domain = linspace(d_beg,d_end,zd);
        target = f(domain);
        dz = [];
        for z_idx = 1 : zd
            zm = rand(xd,yd) < target(z_idx);
            dz(:,:,z_idx) = zm;
        end
        
        if     isequal(dir, [0 0 0])
            dm = rand(xd,yd,zd);
        elseif isequal(dir, [0 0 1])
            dm = dz;
        elseif isequal(dir, [0 1 0])
            dm = dy;
        elseif isequal(dir, [0 1 1])
            dm = dy & dz;
        elseif isequal(dir, [1 0 0])
            dm = dx;
        elseif isequal(dir, [1 0 1])
            dm = dx & dz;
        elseif isequal(dir, [1 1 0])
            dm = dx & dy;
        elseif isequal(dir, [1 1 1])
            dm = dx & dy & dz;
        end
    end

    %% initialize cell types
    function [] = initialize_occupation()
        initial_occupation_pattern;
        occupation = permute(occupation, [2 1 3]);
    end  % function initialize_occupation

    %% initial occupation pattern
    function [] = initial_occupation_pattern()
        occupation(:, 1,             mid_z) = 1;  % place a layer of vessel at the bottom
        occupation(:, 2:mid_y,       mid_z) = 4;  % place a stack of layers of betas just above the vessels
        occupation(:, mid_y+1:y_dim, mid_z) = 3;  % place a stack of layers of alphas just above the alphas
    end  % function initial_occupation_pattern

end  % function simulate_metabolic_symbiosis_2d
