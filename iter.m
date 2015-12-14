clear; 
clc;

% Parameters
N = 10000; % Number of points
n = ceil( N * 0.75 ); % Number of observed uniform points
iflag = 1;
eps = 1E-14;

% Number of iterations
iterations = 1000000;

% The signal
f = @(x) sin( x )+ 2 * cos( 3.*x ) - 3 *cos( 5 .* x );  

% Band limit
c = 7;

% The low pass filter
filter = [  ones( c+1, 1 );
            zeros( N - (2*c+1), 1 );
            ones( c, 1 )  
            ];

% Create a uniform mesh on [-pi, pi) 
x_uni = linspace( -pi, pi ,N + 1 )';
x_uni = x_uni( 1 : end-1);

% Generate a (sorted) random mesh
x_nu = sort( unifrnd( -pi, x_uni(end), [N - 2, 1] ) );
x_nu = [ x_uni(1) ; x_nu ; x_uni(end) ];

% The part of the meshes on which we observe the signal    
x_uni_in = x_uni( 1 : n );
x_nu_in  = x_nu( 1 : n );

% Part of mesh we don't observe
x_uni_out = x_uni( n+1 : end );     
x_nu_out  = x_nu ( n+1 : end );     

% Initialize uniforms
f_uni_in  = f( x_uni_in );                 % Uniform observed mesh
f_uni_out = zeros( length(x_uni_out), 1 ); % Uniform unobserved mesh
f_uni_rec = [ f_uni_in ; f_uni_out ];      % All the uniform mesh

% Initialize non - uniforms
f_nu_in  = f( x_nu_in ); % Non - uniform observed mes
f_nu_rec = [ f_nu_in ; zeros( length(x_nu_out), 1 ) ]; % All the non uniform mesh

% Papoulis Gerchberg iterations.
for i = 1:iterations
    
    %%%%%%%%%%%%% Uniform mesh PG iteration %%%%%%%%%%%%%%

    % Take standard DFT
    f_uni_hat = fft( f_uni_rec );
     
    % Apply low-pass filter to kill high fequencies 
    f_uni_hat = f_uni_hat .* filter;    

    % Use iDFT to reconstruct signal in time domain
    f_uni_rec = ifft( f_uni_hat );
    f_uni_out = f_uni_rec( length( f_uni_in ) + 1:end );
    f_uni_rec = [ f_uni_in ; f_uni_out ];


    %%%%%%%%%%% Non - uniform PG iteration %%%%%%%%%%%%%%%%
    
    % Get only the required bands using type-1 NUFFT. Like using a low-pass filter
    f_nu_hat              = nufft1d1( N, x_nu, f_nu_rec, iflag, eps, 2*c+1 );
    
    % Reconstruct using type-2 NUFFT.
    f_nu_rec( n+1 : end ) = nufft1d2( N - n, x_nu_out,  -iflag, eps, length(f_nu_hat), f_nu_hat );
    
    % When to stop the process because it diverges
    if any(isnan(f_nu_rec)) 
        break
    end
end

% Boring Plotting
figure
plot( x_nu , f_nu_rec+0.01, x_uni, f_uni_rec-0.01, x_uni, f(x_uni) )
legend( {'NUFFT', 'FFT', 'Truth'}, 'Position', [0.15,0.2,0.25,0.1]  ) 
a = 'Reconstruction of ';
b = func2str(f);
cc= ' using ';
d = num2str(i);
e = ' iterations.';
tit = strcat(a,b, cc,d,e);
title( tit );  