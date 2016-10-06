function [] = simulate( simulation )

global simulation_name;
simulation_name = simulation;
temp_path = path;
finish_val = onCleanup( @() post_simulation( temp_path ) );
addpath( 'external' );
addpath( 'util' );
addpath( simulation );
simulator;

    function post_simulation( p )
        path( p );
    end

end  % end function simulate