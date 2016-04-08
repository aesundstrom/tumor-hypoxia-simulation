% global parameters
s_dim               = 40;     % symmetric dimension (used for convenience, when x, y, z dims are identical)
x_dim               = s_dim;  % x dimension
y_dim               = s_dim;  % y dimension
z_dim               = 1;      % z dimension
num_cell_types      = 6;      % number of cell types
num_particle_types  = 4;      % number of particle types
num_iters           = 1000;   % number of simulation clock ticks
d_tau               = 1;      % time scale separation paremeter
plot_every          = 1;      % number of clock ticks between plottings
reproduce_every     = 1;      % number of clock ticks between probablilistic reproductions
delay_occupation_by = 0;      % number of clock ticks before initializing cell type occupations (used to establish gradients prior to exposing cells)
plot_3d             = 1;      % predicate for plotting 3D occupation per cell type

% output predicates
output_results                = 1;          % predicate for outputting plots (global switch for those below)
output_cell_types             = 1;          % predicate for outputting cell types
output_time_series            = 1;          % predicate for outputting time series
output_local_fitness_3d       = 1;          % predicate for outputting local fitness in 3D
output_local_fitness          = 1;          % predicate for outputting local fitness
output_neighborhood_fitness   = 0;          % predicate for outputting neighborhood fitness
output_particle_concentration = 0;          % predicate for outputting particle concentration
output_only_cell_types        = [5];        % subset of cell types to output
output_only_particle_types    = [1 2 3 4];  % subset of particle types to output

% diffusion rate of each particle type
diffusion_rate      = [10   0.1  0.1  1];

% initial concentration of each particle type
init_concentration  = [0.1  0.1  0.1  0];

% basal lower bound of each particle type for each cell type
basal_lower         = [0    0    0    0;   % vessel
                       0    0    0    0;   % empty
                       0    0    0    0;   % alpha
                       0    0    0    0;   % beta
                       0    0    0    0;   % gamma
                       0    0    0    0];  % delta

% basal upper bound of each particle type for each cell type                   
basal_upper         = [Inf  Inf  Inf  Inf;   % vessel
                       Inf  Inf  Inf  Inf;   % empty
                       Inf  Inf  Inf  Inf;   % alpha
                       Inf  Inf  Inf  Inf;   % beta
                       Inf  Inf  Inf  Inf;   % gamma
                       Inf  Inf  Inf  Inf];  % delta
                       
% consume rate of each particle type for each cell type
consume_rate        = [0    0    0    0;   % vessel
                       0    0    0    0;   % empty
                       0    0    0    0;   % alpha
                       0    0    0    0;   % beta
                       0    0    0    0;   % gamma
                       0    0    0    0];  % delta

% release rate of each particle type for each cell type                   
release_rate        = [0    0    0    0;   % vessel
                       0    0    0    0;   % empty
                       0    0    0    0;   % alpha
                       0    0    0    0;   % beta
                       0.2  0.2  0    0;   % gamma
                       0    0    0    0];  % delta
                       
% impact factor of each particle type for each cell type                   
impact_factor       = [0    0    0    0;   % vessel
                       0    0    0    0;   % empty
                       0    0    0    0;   % alpha
                       0    0    0    0;   % beta
                       0    0    0    0;   % gamma
                       0    0    0    0];  % delta

% replacement status of each cell type                   
replaceable         = [0;   % vessel
                       1;   % empty
                       0;   % alpha
                       0;   % beta
                       0;   % gamma
                       0];  % delta

% reproductive status of each cell type                   
reproductive        = [0;   % vessel
                       1;   % empty
                       0;   % alpha
                       0;   % beta
                       1;   % gamma
                       0];  % delta
                   
% color of each cell type                   
cell_type_color_map = [1    1    1;      % vessel (white)
                       0    0    1;      % empty  (blue)
                       1    0    0;      % alpha  (red)
                       0    1    0;      % beta   (green)
                       1    0.5  0;      % gamma  (orange)
                       0.73 0.33 0.83];  % delta  (purple)

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
conditionals(4) = {
                    { { {2 '>' 0.11} } { {'c' 2 0.01} {'r' 3 0.2} } }
                  };  % beta
conditionals(5) = {
                    { { {1 '>' 0.11} } { {'c' 1 0.01} {'i' 1 0.1} } }
                    { { {3 '>' 0.11} } { {'c' 3 0.01} {'i' 3 10} } }
                  };  % gamma
