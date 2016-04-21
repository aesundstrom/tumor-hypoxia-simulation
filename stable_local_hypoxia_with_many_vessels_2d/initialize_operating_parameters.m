% global parameters
global s_dim x_dim y_dim z_dim num_cell_types num_particle_types num_iters d_tau plot_every reproduce_every delay_occupation_by plot_3d;
s_dim               = 40;     % symmetric dimension (used for convenience, when x, y, z dims are identical)
x_dim               = s_dim;  % x dimension
y_dim               = s_dim;  % y dimension
z_dim               = 1;      % z dimension
num_cell_types      = 5;      % number of cell types
num_particle_types  = 1;      % number of particle types
num_iters           = 1000;   % number of simulation clock ticks
d_tau               = 1;      % time scale separation parameter
plot_every          = 1;      % number of clock ticks between plots
reproduce_every     = 1;      % number of clock ticks between probabilistic reproductions
delay_occupation_by = 0;      % number of clock ticks before initializing cell type occupations (used to establish gradients prior to exposing cells)
plot_3d             = 1;      % predicate for plotting 3D occupation per cell type

% output predicates
global output_results output_cell_types output_time_series output_local_fitness_3d output_local_fitness output_neighborhood_fitness output_particle_concentration output_only_cell_types output_only_particle_types;
output_results                = 1;        % predicate for outputting plots (global switch for those below)
output_cell_types             = 1;        % predicate for outputting cell types
output_time_series            = 1;        % predicate for outputting time series
output_local_fitness_3d       = 1;        % predicate for outputting local fitness in 3D
output_local_fitness          = 1;        % predicate for outputting local fitness
output_neighborhood_fitness   = 0;        % predicate for outputting neighborhood fitness
output_particle_concentration = 0;        % predicate for outputting particle concentration
output_only_cell_types        = [3 4];    % subset of cell types to output
output_only_particle_types    = [1];      % subset of particle types to output

% diffusion rate of each particle type
global diffusion_rate;
diffusion_rate      = [0.12];

% initial concentration of each particle type
global init_concentration;
init_concentration  = [0.1];

% basal lower bound of each particle type for each cell type
global basal_lower;
basal_lower         = [0;   % vessel
                       0;   % empty
                       0;   % alpha
                       0;   % beta
                       0];  % gamma

% basal upper bound of each particle type for each cell type                   
global basal_upper;
basal_upper         = [Inf;   % vessel
                       Inf;   % empty
                       Inf;   % alpha
                       Inf;   % beta
                       Inf];  % gamma
                       
% consume rate of each particle type for each cell type
global consume_rate;
consume_rate        = [0;      % vessel
                       0;      % empty
                       0.01;   % alpha
                       0.01;   % beta
                       0];     % gamma

% release rate of each particle type for each cell type                   
global release_rate;
release_rate        = [0.2; % vessel
                       0;   % empty
                       0;   % alpha
                       0;   % beta
                       0];  % gamma
                       
% impact factor of each particle type for each cell type                   
global impact_factor;
impact_factor       = [0;   % vessel
                       0;   % empty
                       1;   % alpha
                       1;   % beta
                       0];  % gamma

% replacement status of each cell type                   
global replaceable;
replaceable         = [0;   % vessel
                       1;   % empty
                       0;   % alpha
                       0;   % beta
                       0];  % gamma

% reproductive status of each cell type                   
global reproductive;
reproductive        = [0;   % vessel
                       1;   % empty
                       1;   % alpha
                       1;   % beta
                       0];  % gamma
                   
% color of each cell type                   
global cell_type_color_map;
cell_type_color_map = [1    1    1;   % vessel (white)
                       0    0    1;   % empty  (blue)
                       1    0    0;   % alpha  (red)
                       0    1    0;   % beta   (green)
                       1    0.5  0];  % gamma  (orange)

% declare conditional triggers and actions:
%     trigger = { <particle type> <comparative operator = {<,=,>}> <value> }
%     action  = { <operator = {a,j,-,+,c,r,i}> <operand 1 = {jump cell type, boolean target value, particle type}> <operand 2 = {particle concentration, impact factor}> }
% note:
%     operator:  number of operands:
%     {a}        0
%     {j,-,+}    1
%     {c,r,i}    2
% example conditional triggers and actions:
%     conditionals(n) = { { { {t1} {t2} {t3} } { {a1} {a2}      } }
%                         { { {t1} {t2}      } { {a1} {a2} {a3} } }
%                       };
global conditionals;
conditionals    = containers.Map('KeyType', 'int32', 'ValueType', 'any');
conditionals(3) = {
                    { { {1 '<' 0.07} } { {'j' 4} } }
                  };  % alpha
conditionals(4) = {
                    { { {1 '<' 0.05} } { {'j' 5} } }
                    { { {1 '>' 0.07} } { {'j' 3} } }
                  };  % beta
