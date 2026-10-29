%% Model parameters

% -----------------------------------
% Numerical parameters
% -----------------------------------
n = 100;    % grid size
h = 1/n;    % spatial step size
k = h;    % time step size

% -----------------------------------
% Physical  parameters
% -----------------------------------
k_0 = .28; %2.136 * 1e-1;     % thermal diffusivity constant
theta = -3*pi/4;                       % wind direction
m = n/8;                            % wind magnitude
v = [m*cos(theta),m*sin(theta)];    % wind vector
Arr = 1.8793 * 1e2;     % empirical burn rate
Bc = 5.5849 * 1e2;      % empircal proportionality coefficient
C = 4.8372 * 1e-5;      % empirical convenction cooling coefficient
% Cs = 1.625 * 1e-1;      % empirical fuel relative disapearrance rate
Cs = 0.5;             % faster fuel relative disapperance rate
Ta = 300;               % ignition temp
IC = 1200;              % initial temp bump

T = 400; % Final time step

%% Construction of spatial discretization matrices

% -----------------------------------
% Laplacian matrix
% -----------------------------------
main = -4*ones(n,1);
off = ones(n,1);

A = spdiags([off main off], -1:1, n, n);

M = [];
for i = 1:n
    M = blkdiag(M,A);
end

[l1, l2] = size(M);

I = spdiags([ones(l1,1) ones(l1,1)],[-n, n], l1,l2);
D = M+I;
D = k_0/(h^2)*D;

% -----------------------------------
% v dot gradient of U 
% -----------------------------------
main2 = zeros(n,1);
off = ones(n,1);

B = spdiags([-1*v(1)*off main2 v(1)*off], -1:1, n, n);

M2 = [];
for i = 1:n
    M2 = sparse(blkdiag(M2,B));
end

I2 = spdiags([-1*v(2)*ones(l1,1) v(2)*ones(l1,1)],[-n, n], l1,l2);
D2 = M2+I2;
D2 = 1/(2*h)*D2;


%% Construction of inverse matrix for implicit time-stepping method

% -----------------------------------
% Backward Euler
% -----------------------------------
l_d = length(D);
Eye = sparse(1:l_d,1:l_d,ones(1,l_d),l_d,l_d);
I_A = (Eye-k*(D-D2));
I_A = sparse(I_A);

% -----------------------------------
% Prefactor the matrix (for speed)
% -----------------------------------
%  
[L_A,U_A,P] = lu(I_A);

%% Construction of initial conditions

c = IC;                     % height of initial heated region
x = linspace(0,1,n);
y = x;
[X,Y] = meshgrid(x,y);

% random spot fires
locations = [.9,.9,.05;.75,.3,.05;.9,.1,.05;...
    .1,.9,.05;.5,.9,.05;.9,.5,.05];  % location and size of random spots
int = FireInitializer(locations,n,IC); % create initial condition

size(int)
% int = IC*(max(X,Y) < 0.2);  % initial heated square region

U_int = reshape(int',n^2,1);


% Initialize U, U_star, U_2_star
U = U_int;
U_star = U;
U_2_star = NaN(size(U_star));

% Define initial fuel distribution
S_int = ones(size(U_int));

% Add randomness to fuel distribution
Ran = rand(size(S_int));
S_int = S_int - (.7*rand(size(S_int)));
S_int = S_int.*(Ran > 0.7) + 0.1*(Ran < 0.7);

% Add fuel break to fuel distribution
% r = n/10; % number of rows to zeros
% S_int(n^2/5:n*(n/5+r)) = zeros(size(S_int(n*10:n*(10+r))));

S = S_int;
S_plot_int = (reshape(S,n,n))';
S_plot = S_plot_int;


U_plot = int;
N = length(U);

F = @(u,S) Arr*(S.*exp(-Bc./(u-Ta)) - C*(u-Ta));

%% Time Stepping the surface

% --------------------------------------
% Backward euler and Godunov splitting
% --------------------------------------

% Plotting/movie setup
xi = linspace(0,1,n);
[X,Y] = meshgrid(xi,xi);
Mov = struct('cdata',[],'colormap',[]);
Vid = VideoWriter('revamped_fire.avi');
open(Vid);

for t = 1:T
    % Update U_star -- reaction part
    U_star = U_star + F(U_star,S).*(U > Ta);
    U_star = sparse(U_star);
    
    % Update U_2_star -- diffusion and advection part
         % U_2_star = I_A\U_star;
    U_2_star = U_A\(L_A\(P*U_star));

    % Update U
    U = U_2_star;
    
    % Update S -- fuel consumption
    S = (S./(1 + k*Cs*exp(-Bc./(U-Ta)))).*(U > Ta) + ...
        S.*(U <= Ta);
    
    % Save U for plotting
    U_plot = (reshape(U,n,n))';
    
    % Save S for plotting
    S_plot = (reshape(S,n,n))';
    
    Z = U_plot;
    contour(X,Y,Z,'fill','on');
    axis([0,1,0,1]);
    caxis([0,4500]);
    colorbar;
    colormap('hot');
	Mov(t) = getframe();
     writeVideo(Vid,Mov(t));
end

close(Vid);

%%

figure
%subplot(1,2,1)
contour(X,Y,S_plot_int,'fill','on');
caxis([0,1])
colorbar
title('Initial fuel distibution')

figure
%subplot(1,2,2)
contour(X,Y,S_plot,'fill','on');
caxis([0,1])
colorbar
title('Final fuel distibution')

% v1 = v(1)*(X./X);
% v2 = v(2)*(Y./Y);
% figure
% quiver(X,Y,v1,v2);
% axis([0,1,0,1])
% title('Wind Field')