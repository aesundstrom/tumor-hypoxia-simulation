function [] = simulate( simulation )

temp_path = path;
addpath( 'external' );
addpath( 'util' );
addpath( simulation );
simulator;
path( temp_path );

end  % end function simulate