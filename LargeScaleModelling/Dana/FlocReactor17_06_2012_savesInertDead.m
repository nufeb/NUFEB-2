% simple program cell growth and shoving in a perfectly mixed reactor
tic
%
clc
close all
clear all
clear java
%
% profile on
%
% PathRoot = [uigetdir '\'];
% PathRoot = 'D:\VLavric\FloconImproved\';
PathRoot = 'D:\FlocVersions\FloconWR\';
% PathRoot = 'D:\FlocWR\';
% PathRoot = 'D:\Cristi\People\Dana_Ofiteru\Model\';
datenow = datestr(now,'dd-mmm-yyyy HH-MM');     % saves each run at the time was created
PathOutput = [PathRoot datenow];
%
mkdir(PathOutput);
%
PathInput = [PathRoot 'SavedStructures'];
PathFcns = PathOutput;
%
javaaddpath('.\java_modif\shoving\build\classes\');   % path to the Java routines
QSobj = shoving.QS();                                 % create a Java object
%
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
% rand('twister', sum(100*clock));
%
nGrowth = 1000;                 % number of growth iterations
nGroupAttach =2;               % frequency of attachment events (number of time steps)
nSingleAttach = 1;              % frequency of attachment events (number of time steps)
DivisionFraction = 0.5;         %
AttachementProbability = 0.25;  % not in use for now
nGroupDetach = 2*nGroupAttach+1; % how often detachment is with respect to attachment
pMaxCellsDetach = 0.25; % fraction of the total floc possible to be detached once
ShowImage = 5; % how often the plots will be saved
%
Nx = 450;       % grid nodes in x direction
Ny = 450;       % grid nodes in y direction
%dx = 1e-6;     % m,
dx = 1.1e-6;
%dx = 1.15e-6;   % m, grid size, compromise  since all the cells have this size (het with EPS 1.34e-6)
%dy = 1e-6;      % m, grid size
dy = 1.1e-6;
%dy = 1.15e-6;   % m, grid size
Lx = Nx*dx;  % m, space length in x; defining a grid/matrix for placing the bacteria before sending to COMSOL
Ly = Ny*dy;  % m, space height in y
%
% properties of each initial bacterium (this can be randomize - either
% position/ no of cells/ type of cells etc)
bac_x = [Lx/2 Lx/2-2*dx  Lx/2+dx  Lx/2+dx Lx/2-dx]; % placing the cells around the center and keeping the positive values for the coordinates
bac_y = [Ly/2 Ly/2+dy    Ly/2-dy  Ly/2+dy Ly/2-dy]; %
%
mass_het = 1e-16;
mass_aob = 1e-16;
mass_nob = 1e-16;
mass_eps = 2.6e-17;   % kg for grid 1e-6 % for having eps density 25 kg/m3
mass_inert = 1.1e-16; % kg for grid 1e-6 % for having eps density 105 kg/m3
%
bac_m   = [mass_het  mass_aob  mass_nob  0          0];          % kg, vector of mass of active cell
bac_e_d = [mass_eps  0         0         mass_eps   mass_inert];   % kg, vector of mass of EPS (associated with het) & solid debries
bac_mshov = bac_m + bac_e_d;                       % kg, particles mass for shoving
%
%       HET    AOB     NOB     EPS     DEAD
bac_s = [1      2       3       4       5  ]; % tag for each floc component
% species clone; keeps track of the descendents of the same bacterium, so all be in a
% cluster (aob & nob only)
bac_c = [0      1       2       0       0  ];
% species properties
bac_mmax = [2*mass_het, 2*mass_aob, 2*mass_nob];      % critical mass before division
%
bac_h    = dx;             % m, cell "height" for 2d approximation
%
bac_rho  = 100;            % kg/m3, biomass cell density
bac_rhod = 105;            % kg/m3, solid debries density (after a cell dies)
bac_rhoe = 25;             % kg/m3, biomass EPS density
%
bac_n    = length(bac_m);  % number of cells
%
dz       = bac_h;          % m
Lbl      = 20*dx;          % m, boundary layer thickness
% Recalculate the maximum colony radius
R_floc   = max(sqrt((bac_x-Lx/2).^2 + (bac_y-Ly/2).^2)) + 10e-6;  % maximum distance + 10 microns
Rcol = [2.5, 2];
%
% Reactor initialization - how much the concentrations are when the process starts
% with innoculum
%
c_o2_reactor_init  = 5e-3;   % kg/m^3
% current o2 concentration in the reactor
c_o2_reactor  = c_o2_reactor_init;
c_o2_sat      = 9e-3;        % kg/m^3
%
c_nh4_reactor_init = 40e-3;  % kg/m^3
% current nh4 concentration in the reactor
c_nh4_reactor = c_nh4_reactor_init;
%
c_no2_reactor_init = 1e-4;   % kg/m^3
% current no2 concentration in the reactor
c_no2_reactor = c_no2_reactor_init;
%
c_s_reactor_init   = 80e-3;  % kg/m^3
% current s concentration in the reactor
c_s_reactor   = c_s_reactor_init;
%
c_no3_reactor_init = 1e-4;   % kg/m^3
% current no3 concentration in the reactor
c_no3_reactor = c_no3_reactor_init;
%
% Reactor's inlet
%
c_o2_reactor_inlet  = c_o2_reactor_init;
c_o2_R  = c_o2_reactor;      % this is the vector of results in time
%
c_nh4_reactor_inlet = c_nh4_reactor_init;
c_nh4_R = c_nh4_reactor;    % this is the vector of results in time
%
c_no2_reactor_inlet = c_no2_reactor_init;
c_no2_R = c_no2_reactor;    % this is the vector of results in time
%
c_s_reactor_inlet   = c_s_reactor_init;
c_s_R   = c_s_reactor;      % this is the vector of results in time
%
c_no3_reactor_inlet = c_no3_reactor_init;
c_no3_R = c_no3_reactor;    % this is the vector of results in time
%
V_reactor = 1;              % m^3
mBiomass  = 0.06;           % kg - innoculum mass
%
% operating parameters of the system
%
alpha = 0.2;    % recycling fraction
beta  = 0.01;   % purje (wastege line) fraction
%
HRT = 7/24;                              % days
SRT = HRT * (alpha+beta)/beta/(1+alpha); % days
%
q_o2 = 0.5;                               % kg / day / m^3 specific flow of oxygen
%
% initial distribution of het:aob:nob is 1:1:1
%
c_het_reactor_init = mBiomass / V_reactor * 1/3;
c_aob_reactor_init = mBiomass / V_reactor * 1/3;
c_nob_reactor_init = mBiomass / V_reactor * 1/3;
%
% diffusion coeficients in water
D_o2  = 2e-9*3600*24;     % m^2/d
D_nh4 = 1.4e-9*3600*24;   % m^2/d
D_no2 = 1.2e-9*3600*24;   % m^2/d
D_s   = 0.5e-9*3600*24;   % m^2/d
D_no3 = 1.2e-9*3600*24;   % m^2/d
% kinetic parameters
mu_max_het = 6;      % 1/d miu_max HET  (Henze ASM2d and others - but consider that here there is no storage
mu_max_aob = 0.76;   % 1/d miu_max AOB (Rittmann and McCarty, 2001) at 20 degrees
mu_max_nob = 0.81;    % 1/d miu_max NOB (Rittmann and McCarty, 2001) at 20 degrees
K_o2_het  = 0.81e-3; % kg/m3
K_s_het   = 10e-3;   % kg/m3
K_no2_het = 0.3e-3;  % kg/m3 % Alpkvist et al., 2006
K_no3_het = 0.3e-3;  % kg/m3
K_o2_aob  = 0.5e-3;  % kg/m3 (Rittmann and McCarty, 2001) at 20 degrees
K_nh4_aob = 1e-3;    % kg/m3 (Rittmann and McCarty, 2001) at 20 degrees
K_o2_nob  = 0.68e-3; % kg/m3 (Rittmann and McCarty, 2001)
K_no2_nob = 1.3e-3;  % kg/m3 (Rittmann and McCarty, 2001) at 20 degrees
% yields
Y_het = 0.61;       % Ni et al., 2009 (Water Research) gCODx/gCODs
Y_eps = 0.18;       % Ni et al., 2009 (Water Research) gCOD_eps/gCODs
Y_eps_het = Y_eps/Y_het; %
Y_s_het   = 1/Y_het;
Y_o2_het  = (1-Y_het-Y_eps)/Y_het;
Y_no2_het = (1-Y_het-Y_eps)/1.71/Y_het;
Y_no3_het = (1-Y_het-Y_eps)/2.86/Y_het;
eta_het   = 0.6;    % reduction factor in anoxic conditions ASM 2d
Y_aob     = 0.33;   %(Rittmann and McCarty, 2001) at 20 degrees
Y_o2_aob  = (3.42-Y_aob)/Y_aob;
Y_nh4_aob = 1/Y_aob;
Y_no2_aob = 1/Y_aob;
Y_nob     = 0.083;   %(Rittmann and McCarty, 2001) at 20 degrees
Y_o2_nob  = (1.15-Y_nob)/Y_nob;
Y_no2_nob = 1/Y_nob;
Y_no3_nob = 1/Y_nob;
%
% decay constants
%
Y_I   = 0.4;  % % Alpkvist et al., 2006
b_het = 0.15; % 1/d Rittman and McCarty, 2001 at 20 degrees
b_nob = 0.11; % 1/d Rittmann and McCarty, 2001 at 20 degrees
b_aob = 0.11; % 1/d Rittmann and McCarty, 2001 at 20 degrees
b_eps = 0.17; % 1/d Ni et al., 2009 (Water Research)
%
% the decay is accounted for at the level of the generic floc, throughout the virtual
% mass of the dead cells which will increase till bigger than the actual mass of het,
% aob or nob cells, moemnt at which death will manifest
%
m_decay_het = 0; % mass of virtual het dead cells
m_decay_aob = 0; % mass of virtual aob dead cells
m_decay_nob = 0; % mass of virtual nob dead cells
%
% a vector in which the number of flocs which will be memorized every step
% here is the first number, as resulted from the innoculum
nFlocs(1) = mBiomass / sum(bac_m+bac_e_d);
%
col         = {[0 0.7 0 ], [0 0 1 ], [1 0 0], [0.7 0.7 0.7], [0 0 0]}; % RGB
k           = 1.1;                                            % factor for overlapping (k>1 for non-overlapping)
%
% the radii of active / EPS and solid debries
% for now, dividing to EPS density only is good, since there is no solid debries
% consecutive to a dead cell but, in the main loop, it will be changed
bac_r       = sqrt((bac_m/bac_rho+bac_e_d/bac_rhoe)/pi/bac_h);  % radius of the cell with eps
bac_ra      = sqrt((bac_m/bac_rho)/pi/bac_h);                   % radius of the active biomass
% attachment parameters
search_dist = 3 * bac_r(1);                                   % search distance around one cell
allowed_overlap = 0.01 * bac_r(1);                            % allowed overlaping before stopping the shoving
R_col = 2.5;                                                  % safety factor for colony size
%
iterShoving  = 5000;
t   = 0;
t_R = t;        % this is the vector storing the integration time
dt  = 0.1;      % day; time step for integration
%
clear fem
% Geometry
g1=rect2(num2str(Lx),num2str(Ly),'base','corner','pos',{'0','0'},'rot','0');
g2=circ2(num2str(R_floc),'base','center','pos',{num2str(Lx/2),num2str(Ly/2)},'rot','0');

% Geometry objects
clear s
s.objs={g1,g2};
s.name={'R1','C1'};
s.tags={'g1','g2'};
fem.draw=struct('s',s);
clear s
fem.geom=geomcsg(fem);
%
% Definition of solutes mass balances with diffusion and reaction
appl.mode.class = 'Diffusion';
appl.dim        = {'c_o2','c_nh4','c_no2','c_s', 'c_no3'};
appl.module     = 'CHEM';
appl.gporder    = 4;
appl.cporder    = 2;
appl.assignsuffix = '_chdi';
%
bnd.c0 = {{'cb_o2';'cb_nh4';'cb_no2';'cb_s';'cb_no3'},0};
bnd.type = {{'C';'C';'C';'C';'C'},{'cont';'cont';'cont';'cont';'cont'}};
bnd.ind = [1,1,1,1,2,2,2,2];
appl.bnd = bnd;
clear bnd
%
equ.D           = {{'D_eff_o2'; 'D_eff_nh4'; 'D_eff_no2'; 'D_eff_s'; 'D_eff_no3'}};
equ.R           = {{'r_o2'; 'r_nh4'; 'r_no2'; 'r_s'; 'r_no3'}};
equ.init        = {{'cb_o2';'cb_nh4';'cb_no2';'cb_s';'cb_no3'}};
equ.ind         = [1,1];
appl.equ        = equ;
clear equ
%
fem.appl{1}     = appl;
clear appl
%
fem.frame       = {'ref'};
fem.border      = 1;
fem.outform     = 'general';
%
units.basesystem= 'SI';
fem.units       = units;
clear units;
%
% Coupling variable elements
clear elemcpl
% Integration coupling variables
clear elem
elem.elem = 'elcplscalar';
elem.g = {'1'};
src = cell(1,1);
clear equ
equ.expr = {{'r_o2*dz'},{'r_s*dz'},{'r_no2*dz'},{'r_no3*dz'},{'r_nh4*dz'}};
equ.ipoints = {{'4'},{'4'},{'4'},{'4'},{'4'}};
equ.frame = {{'ref'},{'ref'},{'ref'},{'ref'},{'ref'}};
equ.ind = {{'1','2'}};
src{1} = {{},{},equ};
elem.src = src;
geomdim = cell(1,1);
geomdim{1} = {};
elem.geomdim = geomdim;
elem.var = {'r_o2_int','r_s_int','r_no2_int','r_no3_int','r_nh4_int'};
elem.global = {'1','2','3','4','5'};
elem.maxvars = {};
elemcpl{1} = elem;
fem.elemcpl = elemcpl;

% Global expressions
fem.globalexpr = {... % each grid cell (x,y) has its own microorganism concentration
  'r_het',    'cbio_het(x,y)*mu_max_het*c_s/(K_s_het+c_s)*c_o2/(K_o2_het+c_o2)+cbio_het(x,y)*eta_het*mu_max_het*c_s/(K_s_het+c_s)*c_no2/(K_no2_het+c_no2)*K_o2_het/(K_o2_het+c_o2)+cbio_het(x,y)*eta_het*mu_max_het*c_s/(K_s_het+c_s)*c_no3/(K_no3_het+c_no3)*K_o2_het/(K_o2_het+c_o2)', ...
  'r_aob',    'cbio_aob(x,y)*mu_max_aob*c_nh4/(K_nh4_aob+c_nh4)*c_o2/(K_o2_aob+c_o2)', ...
  'r_nob',    'cbio_nob(x,y)*mu_max_nob*c_no2/(K_no2_nob+c_no2)*c_o2/(K_o2_nob+c_o2)', ...
  'r_het_no2' 'cbio_het(x,y)*eta_het*mu_max_het*c_s/(K_s_het+c_s)*c_no2/(K_no2_het+c_no2)*K_o2_het/(K_o2_het+c_o2)',...
  'r_het_no3' 'cbio_het(x,y)*eta_het*mu_max_het*c_s/(K_s_het+c_s)*c_no3/(K_no3_het+c_no3)*K_o2_het/(K_o2_het+c_o2)',...
  'r_o2',     '-Y_o2_het*r_het-Y_o2_aob*r_aob-Y_o2_nob*r_nob', ...
  'r_nh4',    '-Y_nh4_aob*r_aob', ...
  'r_no2',    'Y_no2_aob*r_aob-Y_no2_nob*r_nob -Y_no2_het*r_het_no2', ...
  'r_s',      '-Y_s_het*r_het',...
  'r_no3',    'Y_no3_nob*r_nob-Y_no3_het*r_het_no3', ...
  'D_eff_o2', 'f_eff(x,y)*D_o2*(1-0.43*cbio_tot^0.92/(11.19+0.27*cbio_tot^0.99))', ...
  'D_eff_nh4','f_eff(x,y)*D_nh4*(1-0.43*cbio_tot^0.92/(11.19+0.27*cbio_tot^0.99))', ...
  'D_eff_no2','f_eff(x,y)*D_no2*(1-0.43*cbio_tot^0.92/(11.19+0.27*cbio_tot^0.99))', ...
  'D_eff_s',  'f_eff(x,y)*D_s*(1-0.43*cbio_tot^0.92/(11.19+0.27*cbio_tot^0.99))', ...
  'D_eff_no3','f_eff(x,y)*D_no3*(1-0.43*cbio_tot^0.92/(11.19+0.27*cbio_tot^0.99))', ...
  'rhs_o2',   'r_o2_int*nFlocs/V_reactor+(c_o2_reactor_inlet-cb_o2)/HRT+q_o2*(1-cb_o2/c_o2_sat)', ...
  'rhs_s',    'r_s_int*nFlocs/V_reactor+(c_s_reactor_inlet-cb_s)/HRT', ...
  'rhs_no2',  'r_no2_int*nFlocs/V_reactor+(c_no2_reactor_inlet-cb_no2)/HRT', ...
  'rhs_no3',  'r_no3_int*nFlocs/V_reactor+(c_no3_reactor_inlet-cb_no3)/HRT', ...
  'rhs_nh4',  'r_nh4_int*nFlocs/V_reactor+(c_nh4_reactor_inlet-cb_nh4)/HRT', ...
  'cbio_tot', 'cbio_het(x,y)+cbio_aob(x,y)+cbio_nob(x,y)+cbio_eps(x,y)+cbio_dead(x,y)'...
  };
%
% (1-0.43*cbio_tot^0.92/(11.19+0.27*cbio_tot^0.99)) as factor for the effective diffusion in biofilm taken from Renslow et al., 2010
%
% ODE Settings for integrating reactor
clear ode
ode.dim={'cb_o2','cb_s','cb_no2','cb_no3','cb_nh4'};
ode.f={'-cb_o2t+rhs_o2',...
  '-cb_st+rhs_s',...
  '-cb_no2t+rhs_no2',...
  '-cb_no3t+rhs_no3',...
  '-cb_nh4t+rhs_nh4'};
ode.init={'cb0_o2','cb0_s','cb0_no2','cb0_no3','cb0_nh4'};
ode.dinit={'0','0','0','0','0'};
units.basesystem = 'SI';
ode.units = units;
fem.ode=ode;
clear ode units
%
xg = dx/2:dx:Lx;                    % x-grid positions for biomass and solutes distribution
yg = dy/2:dy:Ly;                    % y-grid positions for biomass and solutes distribution
[xgrid ygrid] = meshgrid(xg,yg);
grid = [xgrid(:)'; ygrid(:)'];      % generating the grid
%
% Put cell masses on a square grid => biomass concentration
c_bac = zeros(Nx,Ny,5);      % initialize matrix with biomass HET, AOB, NOB, EPS & inert
for i=1:bac_n                % loop on all the cells
  jj = floor(bac_x(i)/dx)+1; % find index of the cell on the grid
  ii = floor(bac_y(i)/dy)+1;
  %
  if bac_s(i) < 4
    c_bac(ii,jj,bac_s(i)) = bac_m(i); % add cell mass bac_s = 1, 2, 3, 4 or 5
  else
    c_bac(ii,jj,bac_s(i)) = bac_e_d(i);  
  end
end
% Transform biomass in biomass concentration and area in cell fraction
c_bac = c_bac/(dx*dy*dz);
% Create a matrix of diffusion coefficients (very large in bulk liquid, and
% normal in the BL and floc)
f_eff = ones(Ny,Nx)*1e3;
%
% save the grid in the form appropriate for Comsol to read it
%
for i=1:5
  save_gridded_data(c_bac(:,:,i), xg, yg, PathOutput, ['c_bac_' num2str(i) '_'], 1);
end
%
save_gridded_data(f_eff(:,:), xg, yg, PathOutput, 'f_eff_', 1);
%
% Define the Matlab figure
% fh=figure('PaperPosition',[0.01 0.0 3.5 3.5],'Position',[ 100 100 500 500 ],'Color','w','PaperPositionMode','auto');
% hax=axes('Position',[ 0.05 0.05 0.95 0.95 ],'Layer','top');
%
% Growth loop (iterGrowth = counter for growth steps)
%
clear fcns
%
% save initial Reactor data
%
nameReactorFile = [PathOutput '\ReactorResults.txt'];
%
file = fopen(nameReactorFile,'w');
fprintf(file,'%% %15s%15s%15s%15s%15s%15s%15s',...
  't_R', 'nFlocs', 'c_s_R', 'c_o2_R', 'c_nh4_R', 'c_no2_R', 'c_no3_R');
fprintf(file,'\n');
fclose(file);
file = fopen(nameReactorFile,'a');
fprintf(file,'%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e',...
  t_R, nFlocs, c_s_R, c_o2_R, c_nh4_R, c_no2_R, c_no3_R);
fprintf(file,'\n');
fclose(file);

pHET  = length(find(bac_s == 1)) / length(bac_s);
pAOB  = length(find(bac_s == 2)) / length(bac_s);
pNOB  = length(find(bac_s == 3)) / length(bac_s);
pEPS  = length(find(bac_s == 4)) / length(bac_s);
pDEAD = length(find(bac_s == 5)) / length(bac_s);
%
nameFlocStatistics = [PathOutput '\FlocStatistics.txt'];
file = fopen(nameFlocStatistics,'w');
fprintf(file,'%% %15s%15s%15s%15s%15s',...
  'pHET', 'pAOB', 'pNOB', 'pEPS', 'pDEAD');
fprintf(file,'\n');
fclose(file);
file = fopen(nameFlocStatistics,'a');
fprintf(file,'%15.5e%15.5e%15.5e%15.5e%15.5e', pHET, pAOB, pNOB, pEPS, pDEAD);
fprintf(file,'\n');
fclose(file);
%
% save initial running data
%
nameRData = [PathOutput '\RunParameters.txt'];
%
file = fopen(nameRData,'w');
fprintf(file,'%% %15s%15s%15s%15s%15s%15s%15s%15s%15s',...
  'nSingleAttach', 'nGroupAttach', 'nGroupDetach', 'dt', 'q_o2', 'alpha', 'beta', 'HRT', 'SRT');
fprintf(file,'\n');
fclose(file);
file = fopen(nameRData,'a');
fprintf(file,'%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e',...
  nSingleAttach, nGroupAttach, nGroupDetach, dt, q_o2, alpha, beta, HRT, SRT);
fprintf(file,'\n');
fclose(file);
% ==================================================================================
%
% here starts the main integration loop with respect to time
%
for iterGrowth=1:nGrowth
  % Functions
  fcns_method = 'linear';
  %   clear fcns
  fcns{1}.type        ='interp';
  fcns{1}.name        ='cbio_het';
  fcns{1}.method      =fcns_method;
  fcns{1}.extmethod   ='const';
  fcns{1}.filename    =[PathFcns '\c_bac_1_' num2str(iterGrowth,'%03d') '.txt']; % HET cells
  fcns{2}.type        ='interp';
  fcns{2}.name        ='cbio_aob';
  fcns{2}.method      =fcns_method;
  fcns{2}.extmethod   ='const';
  fcns{2}.filename    =[PathFcns '\c_bac_2_' num2str(iterGrowth,'%03d') '.txt']; % AOB cells
  fcns{3}.type        ='interp';
  fcns{3}.name        ='cbio_nob';
  fcns{3}.method      =fcns_method;
  fcns{3}.extmethod   ='const';
  fcns{3}.filename    =[PathFcns '\c_bac_3_' num2str(iterGrowth,'%03d') '.txt']; % NOB cells
  fcns{4}.type        ='interp';
  fcns{4}.name        ='cbio_eps';
  fcns{4}.method      =fcns_method;
  fcns{4}.extmethod   ='const';
  fcns{4}.filename    =[PathFcns '\c_bac_4_' num2str(iterGrowth,'%03d') '.txt']; % EPS cells
  fcns{5}.type        ='interp';
  fcns{5}.name        ='cbio_dead';
  fcns{5}.method      =fcns_method;
  fcns{5}.extmethod   ='const';
  fcns{5}.filename    =[PathFcns '\c_bac_5_' num2str(iterGrowth,'%03d') '.txt']; % DEAD cells
  fcns{6}.type        ='interp';
  fcns{6}.name        ='f_eff';
  fcns{6}.method      =fcns_method;
  fcns{6}.extmethod   ='const';
  fcns{6}.filename    =[PathFcns '\f_eff_' num2str(iterGrowth,'%03d') '.txt']; % param for differentiate diffusion
  %
  fem.functions = fcns;
  clear fcns
  %
  fem.const = {...
    'D_o2',      num2str(D_o2), ...
    'D_nh4',     num2str(D_nh4), ...
    'D_no2',     num2str(D_no2), ...
    'D_s',       num2str(D_s), ...
    'D_no3',     num2str(D_no3), ...
    'mu_max_het',num2str(mu_max_het), ...
    'mu_max_aob',num2str(mu_max_aob), ...
    'mu_max_nob',num2str(mu_max_nob), ...
    'K_o2_het',  num2str(K_o2_het), ...
    'K_s_het',   num2str(K_s_het), ...
    'K_no2_het'  num2str(K_no2_het), ...
    'K_no3_het'  num2str(K_no3_het), ...
    'K_o2_aob',  num2str(K_o2_aob), ...
    'K_nh4_aob', num2str(K_nh4_aob), ...
    'K_o2_nob',  num2str(K_o2_nob), ...
    'K_no2_nob', num2str(K_no2_nob), ...
    'Y_het',     num2str(Y_het), ...
    'Y_eps',     num2str(Y_eps), ...
    'Y_s_het',   num2str(Y_s_het), ...
    'Y_o2_het',  num2str(Y_o2_het), ...
    'Y_no2_het', num2str(Y_no2_het), ...
    'Y_no3_het', num2str(Y_no3_het), ...
    'eta_het',   num2str(eta_het), ...
    'Y_aob',     num2str(Y_aob), ...
    'Y_o2_aob',  num2str(Y_o2_aob), ...
    'Y_nh4_aob', num2str(Y_nh4_aob), ...
    'Y_no2_aob', num2str(Y_no2_aob), ...
    'Y_nob',     num2str(Y_nob), ...
    'Y_o2_nob',  num2str(Y_o2_nob), ...
    'Y_no2_nob', num2str(Y_no2_nob),...
    'Y_no3_nob', num2str(Y_no3_nob),...
    'cb0_o2',    num2str(c_o2_reactor), ...
    'cb0_nh4',   num2str(c_nh4_reactor), ...
    'cb0_no2',   num2str(c_no2_reactor), ...
    'cb0_s',     num2str(c_s_reactor),...
    'cb0_no3',   num2str(c_no3_reactor),...
    'nFlocs',    num2str(nFlocs(end)), ...
    'V_reactor', num2str(V_reactor), ...
    'alpha',     num2str(alpha), ...
    'HRT',       num2str(HRT), ...
    'q_o2',      num2str(q_o2), ...
    'c_o2_sat',  num2str(c_o2_sat), ...
    'c_no2_reactor_inlet',  num2str(c_no2_reactor_inlet), ...
    'c_o2_reactor_inlet',   num2str(c_o2_reactor_inlet), ...
    'c_no3_reactor_inlet',  num2str(c_no3_reactor_inlet), ...
    'c_nh4_reactor_inlet',  num2str(c_nh4_reactor_inlet), ...
    'c_s_reactor_inlet',    num2str(c_s_reactor_inlet), ...
    'dz',       num2str(dz)...
    };
  
  % Geometry
  % this is the rectangle where the floc develops
  g1=rect2(num2str(Lx),num2str(Ly),'base','corner','pos',{'0','0'},'rot','0');
  % this is the circle surrounding the floc for a finer meshing
  g2=circ2(num2str(R_floc),'base','center','pos',{num2str(Lx/2),num2str(Ly/2)},'rot','0');
  
  % Geometry objects
  clear s
  s.objs={g1,g2};
  s.name={'R1','C1'};
  s.tags={'g1','g2'};
  fem.draw=struct('s',s);
  clear s
  fem.geom=geomcsg(fem);
  
  %
  % Multiphysics
  fem=multiphysics(fem);
  % Initialize mesh
  fem.mesh=meshinit(fem, ...
    'hauto',9, ...
    'hmaxsub',[2,5e-6]);
  % Extend mesh
  fem.xmesh=meshextend(fem, 'dofversion',1);
  %
  % Solve problem
  if iterGrowth==1
    fem.sol=femtime(fem, ...
      'solcomp',{'c_no2','c_nh4','c_no3','c_s','cb_nh4','c_o2','cb_o2','cb_no2','cb_s','cb_no3'}, ...
      'outcomp',{'c_no2','c_nh4','c_no3','c_s','cb_nh4','c_o2','cb_o2','cb_no2','cb_s','cb_no3'}, ...
      'blocksize','auto', ...
      'odesolver','genalpha', ...
      'tlist',t:dt/4:t+dt, ...
      'tout','tlist', ...
      'incrdelay','off', ...
      'nlsolver','manual', ...
      'ntolfact',1, ...
      'maxiter',4, ...
      'dtech','const', ...
      'damp',1.0, ...
      'jtech','minimal');
  else
    sol = asseminit(fem,'init',fem0);
    fem.sol=femtime(fem, ...
      'init',sol, ...
      'solcomp',{'c_no2','c_nh4','c_no3','c_s','cb_nh4','c_o2','cb_o2','cb_no2','cb_s','cb_no3'}, ...
      'outcomp',{'c_no2','c_nh4','c_no3','c_s','cb_nh4','c_o2','cb_o2','cb_no2','cb_s','cb_no3'}, ...
      'blocksize','auto', ...
      'odesolver','genalpha', ...
      'tlist',t:dt/4:t+dt, ...
      'tout','tlist', ...
      'incrdelay','off', ...
      'nlsolver','manual', ...
      'ntolfact',1, ...
      'maxiter',4, ...
      'dtech','const', ...
      'damp',1.0, ...
      'jtech','minimal');
  end
  
  fem0 = fem;
  %
  % Save fem structure to be used in the file graphs_2d.m
  %
  save([PathOutput '\fem.' num2str(iterGrowth,'%03d') '.mat'],'fem');
  %
  % Save cells to be used in the file graphs_2d.m
  file = fopen([PathOutput '\cells.' num2str(iterGrowth,'%03d') '.txt'],'wt');
  fprintf(file,'%15s%15s%15s%15s%15s%15s%15s%15s%15s\n',...
    'index','x[um]','y[um]','species','colony','mass_a[kg]','mass_e or d[kg]','rad_eps[um]','rad_act[um]');
  for i=1:bac_n
    fprintf(file,'%15d%15.2f%15.2f%15d%15d%15.5e%15.5e%15.2f%15.2f\n',...
      i, bac_x(i)*1e6, bac_y(i)*1e6, bac_s(i), bac_c(i), bac_m(i),  bac_e_d(i), bac_r(i)*1e6, bac_ra(i)*1e6); % lenght in micrometers
  end
  fclose(file);
  %
  % Put on a grid the variables to be used for biomass growth
  %
  % concentrations for the floc in steady state
  c_s    = reshape(postinterp(fem,'c_s',grid),size(xgrid));   c_s(c_s<0)     = 0;
  c_nh4  = reshape(postinterp(fem,'c_nh4',grid),size(xgrid)); c_nh4(c_nh4<0) = 0;
  c_no2  = reshape(postinterp(fem,'c_no2',grid),size(xgrid)); c_no2(c_no2<0) = 0;
  c_o2   = reshape(postinterp(fem,'c_o2',grid),size(xgrid));  c_o2(c_o2<0)   = 0;
  c_no3  = reshape(postinterp(fem,'c_no3',grid),size(xgrid)); c_no3(c_no3<0) = 0;
  %
  % concentrations in the bioreactor
  cb_v = posteval(fem,'cb_o2','T',t+dt);   cb_o2 = cb_v.d(1);   c_o2(c_o2<0)   = 0;
  cb_v = posteval(fem,'cb_no2','T',t+dt);  cb_no2 = cb_v.d(1);  c_no2(c_no2<0) = 0;
  cb_v = posteval(fem,'cb_no3','T',t+dt);  cb_no3 = cb_v.d(1);  c_no3(c_no3<0) = 0;
  cb_v = posteval(fem,'cb_nh4','T',t+dt);  cb_nh4 = cb_v.d(1);  c_nh4(c_nh4<0) = 0;
  cb_v = posteval(fem,'cb_s','T',t+dt);    cb_s = cb_v.d(1);    c_s(c_s<0)     = 0;
  
  %%% Bacterial attachment
  %
  bac_x = bac_x - Lx/2; % translation of coordinates just for the attachment
  bac_y = bac_y - Ly/2; % translation of coordinates just for the attachment
  %
  [bac_theta,bac_R] = cart2pol(bac_x,bac_y);  % Compute polar coordinates for each cell
  Rmax = max(bac_R);
  Rbl  = Rmax + Lbl;    % boundary layer radius (calculated here but not used in attachment)
  tic
  %
  mFlocOld = sum(bac_m+bac_e_d);
  %
  if mod(iterGrowth,nSingleAttach)==0
    % generate randomly a start position for the new cell
    fi_atta = rand*2*pi;
    R_atta = Rmax + 5*bac_r(1);
    %
    while true
      R_atta = R_atta - bac_r(1); % the attaching cell is moving with bac_r(1) step
      [x_atta,y_atta] = pol2cart(fi_atta,R_atta);
      [dummy_theta,dist] = cart2pol(bac_x - x_atta, bac_y - y_atta); % simlply to speed-up
      if any((dist - 3*bac_r)<0) % reads all the vector bac_r at one time, without looping bt bac_n
        break;
      end
    end
    %
    % create a new cell
    bac_n = bac_n + 1;
    bac_x(bac_n) = x_atta;
    bac_y(bac_n) = y_atta;
    %
    if rand < 0.6
      bac_s(bac_n) = 1;
      bac_c(bac_n) = 0;
      bac_m(bac_n)   = mass_het;
      bac_e_d(bac_n) = mass_eps;
    else           % here would attach only AOB and NOB
      if rand > 0.5
        bac_s(bac_n) = 2; % AOB
        bac_m(bac_n)   = mass_aob;
      else
        bac_s(bac_n) = 3; % NOB
        bac_m(bac_n)   = mass_nob;
      end
      bac_c(bac_n) = max(bac_c)+1;
      bac_e_d(bac_n) = 0;   % AOB and NOB don't have EPS
      
    end
    %
    bac_ra(bac_n)  = sqrt(bac_m(bac_n)/pi/bac_rho/bac_h);
    % here the attached cell is alive, so bac_e_d would have only EPS
    bac_r(bac_n)   = sqrt((bac_m(bac_n)/bac_rho+bac_e_d(bac_n)/bac_rhoe)/pi/bac_h);
    %
  end % attaching a single cell
  %
  %%%%%%%%%%%%%%%%%
  % attach a group of cells
  % load a file with cells attributes
  %
  if mod(iterGrowth,nGroupAttach)==0
    D = dir([PathInput '\*.txt']); % reads the .txt files from the directory
    %
    [a b]= size(D);       % a gives the number of .txt files in the directory; a it's an integer; D it's a structure with 4 fields
    %
    upload = randi(a);% pick randomly one of the .txt files from the directory
    %
    upload_file = D(upload).name;                        % reads the name of the files
    file = fopen([PathInput '\' upload_file],'r');       % open the file
    line = fgets(file);
    %
    [cells,count] = fscanf(file,'%d%f%f%d%d%e%e%f%f',[9,inf]);
    fclose(file);
    %
    [N1 bacatt_n] = size(cells);
    bacatt_x  = cells(2,:)*1e-6;  % m
    bacatt_y  = cells(3,:)*1e-6;  % m
    %
    bacatt_x = bacatt_x - Lx/2;   % puting them back to be centered in (0,0)
    bacatt_y = bacatt_y - Ly/2;
    %
    bacatt_s  = cells(4,:);
    bacatt_c  = cells(5,:);
    bacatt_m  = cells(6,:);
    bacatt_e  = cells(7,:);
    bacatt_r  = cells(8,:)*1e-6;
    bacatt_ra = cells(9,:)*1e-6;
    bacatt_nc = sum(bacatt_c>0);  % number of clones (only for AOB and NOB)
    %
    % center of mass for the coming group of cells
    mtotatt  = sum(bacatt_m);
    mxtotatt = sum(bacatt_m.*bacatt_x);
    mytotatt = sum(bacatt_m.*bacatt_y);
    %
    xcatt = mxtotatt/mtotatt;
    ycatt = mytotatt/mtotatt;
    %Compute polar coordinates for the center of mass of the coming group
    [bacgroup_theta, bacgroup_R] = cart2pol(xcatt,ycatt);
    %
    % generate the random angle for attachement and identifying the cell
    % from the old group to attach to; it is an iterative process until a
    % cell is found in the +-10% region of the generated angle
    %
    ind_theta =[];
    %
    while isempty(ind_theta)
      fi_atta = (1-rand*2)*pi;   % generate a random angle to attach
      %
      % calculates the difference between the angles of the cells in the old
      % group and the attachment angle of the new group
      dif_theta  = bac_theta - fi_atta;
      % selects for the closest cells regarding the angle; now searches in +-10% region from the angle
      ind_theta  = find(abs(dif_theta) <= 0.1*abs(fi_atta));
    end
    % identifies the index of the closest cell from the old group considering the radius
    closest      = find(bac_R == max(bac_R(ind_theta)), 1);
    reffer_theta = bac_theta(closest);
    reffer_R     = bac_R(closest);
    %the carthesian coordinates for the refference cell from the old group,
    % towards which the new group will be moving
    [x_reffer,y_reffer] = pol2cart(reffer_theta,reffer_R);
    % forcing the center of mass for the group to have the same angle as the one for attachment
    bacgroup_theta = reffer_theta;
    %
    % will be given by the radius of the existing floc + radius of the newcommers
    bacgroup_R     = bacgroup_R + Rmax;
    [xcatt_new,ycatt_new] = pol2cart(bacgroup_theta,bacgroup_R);  % recalculating the carthesian coordinates for the center
    %
    dxatt = xcatt_new - xcatt;                             % the distance on 0x to move all the cells;
    dyatt = ycatt_new - ycatt;                             % the distance on 0y to move all the cells;
    %
    bacatt_x = bacatt_x + dxatt;                           % these are the new coordinates for the cells of the coming group
    bacatt_y = bacatt_y + dyatt;
    %
    % checks now how far away are the cells of the
    % coming group from the refference cell
    %
    while true
      %recalculating the carthesian coordinates for the center
      [xcatt_trans,ycatt_trans] = pol2cart(bacgroup_theta,bacgroup_R);
      %
      dxatt = xcatt_trans - xcatt_new;    % the distance on 0x to move all the cells;
      dyatt = ycatt_trans - ycatt_new;    % the distance on 0y to move all the cells;
      %
      xcatt_new = xcatt_trans;
      ycatt_new = ycatt_trans;
      % these are the new coordinates for the cells of the coming group
      bacatt_x = bacatt_x + dxatt;
      bacatt_y = bacatt_y + dyatt;
      % calculates the distance between all the cells
      [dummy_theta,dist_group] = cart2pol(bacatt_x - x_reffer, bacatt_y - y_reffer);
      % of the coming group and the refference cell
      if any((dist_group - 3*bac_r(closest))< 0)
        break;
      end
      % the center of the attaching group is moving with bac_r(closest) step
      bacgroup_R = bacgroup_R - bac_r(closest);
    end
    %the new cells must be added to the old ones
    for q = 1:bacatt_n
      bac_n = bac_n + 1;
      bac_x(bac_n) = bacatt_x(q);
      bac_y(bac_n) = bacatt_y(q);
      bac_s(bac_n) = bacatt_s(q);
      %
      if bac_s(bac_n)==2 || bac_s(bac_n)==3
        bac_c(bac_n) = max(bac_c) + 1;
      else
        bac_c(bac_n) = 0;
      end
      %
      bac_m(bac_n)  = bacatt_m(q);
      bac_e_d(bac_n)  = bacatt_e(q);
      bac_ra(bac_n) = bacatt_ra(q);
      bac_r(bac_n)  = bacatt_r(q);
    end
  end  % attaching a group of cells
  %
  %%%%%%%%%%%%%%%%%
  % detach a group of cells
  %%%%%%%%%%%%%%%%%
  %
  ndetach = 0;
  if mod(iterGrowth,nGroupDetach)==0
    %
    [bac_theta,bac_R] = cart2pol(bac_x,bac_y);  % Compute polar coordinates for each cell
    %
    % find the start angle of detachment
    %
    ind_theta =[];
    %
    while isempty(ind_theta)
      fi_atta = (1-rand*2)*pi;   % generate a random angle to attach
      %
      % calculates the difference between the angles of the cells in the old
      % group and the attachment angle of the new group
      dif_theta  = bac_theta - fi_atta;
      % selects for the closest cells regarding the angle; now searches in +-10% region from the angle
      ind_theta  = find(abs(dif_theta) <= 0.2*abs(fi_atta));
    end
    %
    % find the radius from which the detachment will begin
    % assume that detachment cannot start from the core, but from 75% up of the floc
    R_min_det = max(bac_R(ind_theta)) * (0.5+rand*0.5);
    %
    % find the index of all cells, from the previous cells in the sector, which are
    % in this particular area
    % this keeps only the indexes of the cells which are between the angles and with
    % the radius greater or equal than R_min_det
    %
    % avoid having duplicates in the index vector of removing cells
    %
    ind = unique(ind_theta(bac_R(ind_theta)>=R_min_det)); clear ind_theta
    %
    % find indexes of all AOB & NOB cells from the floc
    %
    indAOBNOBall = sort([find(bac_s == 2) find(bac_s == 3)]);
    %
    % find indexes of aob cells to be detached in the vector of whole cells
    %
    indAOBNOBDetached = intersect(ind, indAOBNOBall);
    %
    % now, add the aob from the collony, if the selected aob belongs to one
    % search if the current aob cell for detachment is in a group
    %
    for iGD = 1:length(indAOBNOBDetached)
      %
      % find if the present aob belongs to a cluster
      %
      ind_cluster = find(bac_c(indAOBNOBDetached(iGD)) == bac_c);
      %
      if length(ind_cluster) > 1
        % guilty of being in a group, so add the whole group indexes
        ind(end+1:end+length(ind_cluster)) = ind_cluster(1:end);
      end
      %
    end
    %
    ind = unique(ind);  ndetach = length(ind);
    %
    bac_m(ind) = [];
    bac_c(ind) = [];
    bac_x(ind) = [];
    bac_y(ind) = [];
    bac_s(ind) = [];
    bac_e_d(ind) = [];
    bac_R(ind) = [];
    Rmax = max(bac_R);
    Rbl  = Rmax + Lbl; % boundary layer radius (calculated here but not used in attachment)
    bac_n = length(bac_m);
    %
    clear ind indAOBNOBDetached indAOBNOBall
  end
  %
  % %%%%%% Detecting and removing islands
  %
  clear bac_R bac_theta
  [bac_theta,bac_R] = cart2pol(bac_x,bac_y);  % Compute polar coordinates for each cell
  %
  bac_x = bac_x + Lx/2;   % going back to the initial coordinates
  bac_y = bac_y + Ly/2;
  %
  ndetisland = 0;
  if mod(iterGrowth,nGroupDetach)==1
    % Mark cells with 1 and liquid with 0
    land = zeros(Ny,Nx);
    for i=1:bac_n                   % loop on all the cells
      jj = floor(bac_x(i)/dx)+1;    % find index of the cell on the grid
      ii = floor(bac_y(i)/dy)+1;
      land(ii,jj)=1;
    end
    %
    %Detecting islands
    %
    % First, mark as land (3) the cells in the middle of the domain
    % (within a square of 6x6 about the center)
    for i=Ny/2-5:Ny/2+5
      for j=Nx/2-5:Nx/2+5
        if land(i,j)==1
          land(i,j)=3;
        end;
      end;
    end;
    %
    % Second, mark as land (3) all the cells in connection with the core
    % (nearest or next-nearest neighbour)
    changes = 1;
    while changes ~= 0
      changes = 0;
      for i=2:Ny-1
        for j=2:Nx-1
          if land(i,j)==1 && ...
              (land(i+1,j)==3 || land(i-1,j)==3 || land(i+2,j)==3 || land(i-2,j)==3 ||...
              land(i,j-1)==3 || land(i,j+1)==3 || land(i,j-2)==3 || land(i,j+2)==3 ||...
              land(i+1,j+1)==3 || land(i-1,j+1)==3 || land(i+2,j+1)==3 || land(i-2,j+1)==3 ||...
              land(i+1,j-1)==3 || land(i-1,j-1)==3 || land(i+2,j-1)==3 || land(i-2,j-1)==3 ||...
              land(i+1,j+2)==3 || land(i-1,j+2)==3 || land(i+2,j+2)==3 || land(i-2,j+2)==3 ||...
              land(i+1,j-2)==3 || land(i-1,j-2)==3 || land(i+2,j-2)==3 || land(i-2,j-2)==3)
            land(i,j)=3;
            changes=changes+1;
          end;
        end;
      end;
    end;
    %
    %Removing islands (all unmarked cells that kept land=1 instead of land=3)
    %
    ind = [];
    for i=1:bac_n                   % loop on all the cells
      jj = floor(bac_x(i)/dx)+1;  % find index of the cell on the grid
      ii = floor(bac_y(i)/dy)+1;
      if land(ii,jj)==1           % if on island, mark the cell for removal
        ind = [ind i];
      end
    end
    %
    if ~isempty(ind)
      bac_m(ind) = [];
      bac_c(ind) = [];
      bac_x(ind) = [];
      bac_y(ind) = [];
      bac_s(ind) = [];
      bac_e_d(ind) = [];
      bac_theta(ind) = [];
      bac_R(ind) = [];
      ndetisland=length(ind);
    end
    %
    Rmax = max(bac_R);
    Rbl  = Rmax + Lbl; % boundary layer radius (calculated here but not used in attachment)
    bac_n = length(bac_m);
    %
    clear land;
  end
  %
  % recalculate the mass of the floc after attachement/detachment
  %
  mFlocNew = sum(bac_m+bac_e_d); % bac_m has the new value of the number of flocs
  %
  % this is how we take into consideration the attachemnet without compromising the existing mass of the cells
  nFlocs(end) = nFlocs(end) * mFlocOld / mFlocNew;
  %
  time_att=toc;
  %
  %%% Growth
  tic
  %
  rx = zeros(1,bac_n); % rate of growth for het, aob & nob
  re = rx; % rate of growth for eps
  rd = rx; % rate of decay for het, aob, nob and inert
  %
  for i=1:bac_n
    % find index of the cell on the grid;
    % they will use the substrates concentration corresponding to the cell grid they belong
    jj = floor(bac_x(i)/dx)+1;
    ii = floor(bac_y(i)/dy)+1;
    %
    switch bac_s(i)
      case 1
        % substrates concentrations are in a matrix now
        rx(i) = bac_m(i)*mu_max_het*c_s(ii,jj)/(K_s_het+c_s(ii,jj))*c_o2(ii,jj)/(K_o2_het+c_o2(ii,jj)) + ...
          bac_m(i)*eta_het*mu_max_het*c_s(ii,jj)/(K_s_het+c_s(ii,jj))*c_no2(ii,jj)/(K_no2_het+c_no2(ii,jj))*K_o2_het/(K_o2_het+c_o2(ii,jj)) + ...
          bac_m(i)*eta_het*mu_max_het*c_s(ii,jj)/(K_s_het+c_s(ii,jj))*c_no3(ii,jj)/(K_no3_het+c_no3(ii,jj))*K_o2_het/(K_o2_het+c_o2(ii,jj)); % kg/time
        re(i) = rx(i)*Y_eps_het; % by default, re will be left zero for solid debries, since bac_s(i) == 5
        rd(i) = b_het * bac_m(i);
      case 2
        rx(i) = bac_m(i)*mu_max_aob*c_nh4(ii,jj)/(K_nh4_aob+c_nh4(ii,jj))*c_o2(ii,jj)/(K_o2_aob+c_o2(ii,jj));
        rd(i) = b_aob * bac_m(i);
      case 3
        rx(i) = bac_m(i)*mu_max_nob*c_no2(ii,jj)/(K_no2_nob+c_no2(ii,jj))*c_o2(ii,jj)/(K_o2_nob+c_o2(ii,jj));
        rd(i) = b_nob * bac_m(i);
      case 4
        rd(i) = b_eps * bac_e_d(i);
    end
  end
  %
  %
  % Find indices of het, aob and nob
  %
  ind_het = find(bac_s == 1); % heterotrophs
  ind_aob = find(bac_s == 2);
  ind_nob = find(bac_s == 3);
  ind_eps = find(bac_s == 4);
  %
  % overall consumptions for one floc at Reactor level
  % (sum of all net consumption rates in one floc)
  %   r_o2_reactor  = -Y_o2_het*sum(rx(ind_het))-Y_o2_aob*sum(rx(ind_aob))-Y_o2_nob*sum(rx(ind_nob));
  %   r_nh4_reactor = -Y_nh4_aob*sum(rx(ind_aob));
  %   r_no2_reactor = Y_no2_aob*sum(rx(ind_aob))-Y_no2_nob*sum(rx(ind_nob))-Y_no2_het*sum(rx(ind_het));
  %   r_s_reactor   = -Y_s_het*sum(rx(ind_het));
  %   r_no3_reactor = Y_no3_nob*sum(rx(ind_nob))-Y_no3_het*sum(rx(ind_het));
  %
  % here, decay is taken care of
  %
  % decay rates for het, aob & nob
  r_decay_het = sum(rd(ind_het));
  r_decay_aob = sum(rd(ind_aob));
  r_decay_nob = sum(rd(ind_nob));
  % mass of the virtual dead cells
  m_decay_het = m_decay_het + r_decay_het * dt;
  m_decay_aob = m_decay_aob + r_decay_aob * dt;
  m_decay_nob = m_decay_nob + r_decay_nob * dt;
  %
  % test if any of masses becomes bigger than the mass of a single corresponding cell
  %
  % total mass of substrate released from the dead cells
  % this will enter the mass balance fr substrate at the bioreactor level
  m_substrate = sum(rd(ind_eps)) * dt;
  bac_e_d(ind_eps) = bac_e_d(ind_eps) - rd(ind_eps)* dt;
  %
  if m_decay_het >= mass_het
    %
    % time for some het cells to die
    %
    while m_decay_het >= mass_het && ~isempty(ind_het)
      % pick-up randomly the het cell to die
      decay_index = randi(length(ind_het));
      % verify if this cell was not already removed
      if bac_m(ind_het(decay_index)) == 0
        continue % this cell is already dead and taken care of
      end
      % eliminate this cell - Y_I will be solid debries, (1-Y_I) will be substrate
      % first, decrease the virtual dead cells mass - one of them becomes real
      m_decay_het = m_decay_het - bac_m(ind_het(decay_index));
      % then, increase the mass of the substrate to be released into the bioreactor
      m_substrate = m_substrate + bac_m(ind_het(decay_index)) * (1-Y_I);
      % the whole eps belonging to this het cell will become substrate
      m_substrate = m_substrate + bac_e_d(ind_het(decay_index));
      % the rest is for the solid debries which remain in the floc taking the place
      % of the late eps
      bac_e_d(ind_het(decay_index)) = bac_m(ind_het(decay_index)) * Y_I;
      % now, remove this het cell from the livings
      bac_m(ind_het(decay_index)) = 0;
      % remove the reaction rate for it
      rx(ind_het(decay_index)) = 0;
      % mark this place as solid debries
      bac_s(ind_het(decay_index)) = 5;
      %
      ind_het(decay_index) = [];
    end
    %
  end
  % now, take care f the aob cells
  if m_decay_aob >= mass_aob
    %
    % time for some aob cells to die
    %
    while m_decay_aob >= mass_aob && ~isempty(ind_aob)
      % pick-up randomly the het cell to die
      decay_index = randi(length(ind_aob));
      % verify if this cell was not already removed
      if bac_m(ind_aob(decay_index)) == 0
        continue % this cell is already dead and taken care of
      end
      % eliminate this cell - Y_I will be solid debries, (1-Y_I) will be substrate
      % first, decrease the virtual dead cells mass - one of them becomes real
      m_decay_aob = m_decay_aob - bac_m(ind_aob(decay_index));
      % then, increase the mass of the substrate to be released into the bioreactor
      m_substrate = m_substrate + bac_m(ind_aob(decay_index)) * (1-Y_I);
      % the rest is for the solid debries which remain in the floc tooking the place
      % of the late eps
      bac_e_d(ind_aob(decay_index)) = bac_m(ind_aob(decay_index)) * Y_I;
      % now, remove this het cell from the livings
      bac_m(ind_aob(decay_index)) = 0;
      % remove the reaction rate for it
      rx(ind_aob(decay_index)) = 0;
      % mark this place as solid debries
      bac_s(ind_aob(decay_index)) = 5;
      %
      ind_aob(decay_index) = [];
    end
    %
  end
  % now, take care f the nob cells
  if m_decay_nob >= mass_nob
    %
    % time for some aob cells to die
    %
    while m_decay_nob >= mass_nob && ~isempty(ind_nob)
      % pick-up randomly the het cell to die
      decay_index = randi(length(ind_nob));
      % verify if this cell was not already removed
      if bac_m(ind_nob(decay_index)) == 0
        continue % this cell is already dead and taken care of
      end
      % eliminate this cell - Y_I will be solid debries, (1-Y_I) will be substrate
      % first, decrease the virtual dead cells mass - one of them becomes real
      m_decay_nob = m_decay_nob - bac_m(ind_nob(decay_index));
      % then, increase the mass of the substrate to be released into the bioreactor
      m_substrate = m_substrate + bac_m(ind_nob(decay_index)) * (1-Y_I);
      % the rest is for the solid debries which remain in the floc tooking the place
      % of the late eps
      bac_e_d(ind_nob(decay_index)) = bac_m(ind_nob(decay_index)) * Y_I;
      % now, remove this het cell from the livings
      bac_m(ind_nob(decay_index)) = 0;
      % remove the reaction rate for it
      rx(ind_nob(decay_index)) = 0;
      % mark this place as solid debries
      bac_s(ind_nob(decay_index)) = 5;
      %
      ind_nob(decay_index) = [];
    end
    %
  end
  %
  %
  % eliminating the EPS granule which is completly disolved
  %
  ind_eps_neg = find(bac_e_d(ind_eps) <= 0); % the indices for "negative" EPS granule
  if ~isempty(ind_eps_neg)
    bac_e_d(ind_eps(ind_eps_neg)) =[];
    bac_s(ind_eps(ind_eps_neg)) = [];
    bac_c(ind_eps(ind_eps_neg)) = [];
    bac_m(ind_eps(ind_eps_neg)) = [];
    bac_x(ind_eps(ind_eps_neg)) = [];
    bac_y(ind_eps(ind_eps_neg)) = [];
    bac_n = length(bac_m);
  end
  %
  % here ends the decay treatment
  %
  clear ind_het ind_aob ind_nob ind_eps
  %
  % balance for the number of flocs which diminishes due to the exit flow
  % here we consider that the pipes and separator have zero volume
  %
  nFlocs(end+1) = nFlocs(end)*(1 - dt/SRT);
  %
  c_o2_R = [c_o2_R; cb_o2];  % stores the concentrations each step
  %
  c_nh4_R = [c_nh4_R; cb_nh4];
  %
  c_no2_R = [c_no2_R; cb_no2];
  %
  c_s_reactor = cb_s + m_substrate * nFlocs(end) / V_reactor;
  c_s_R = [c_s_R; c_s_reactor];
  %
  c_no3_R = [c_no3_R; cb_no3];
  %
  % Biomass mass and size update (Euler)
  bac_m    = bac_m + rx*dt;                % update new biomass
  bac_e_d  = bac_e_d + re*dt;              % update new EPS
  bac_ra   = sqrt(bac_m/pi/bac_rho/bac_h); % the radius of the active biomass
  %
  % compute radii based upon the component occupying the location
  bac_r = zeros(1, bac_n);
  % start with heterotrophs and associated EPS
  ind_eps = [find(bac_s == 1) find(bac_s == 4)];
  % radii of het + EPS
  bac_r(ind_eps) = sqrt((bac_m(ind_eps)/bac_rho+bac_e_d(ind_eps)/bac_rhoe)/pi/bac_h);
  clear ind_eps
  %
  % continue with aob & nob
  ind_aob_nob = [find(bac_s == 2) find(bac_s == 3)];
  % radii of aob & nob
  bac_r(ind_aob_nob) = sqrt(bac_m(ind_aob_nob)/bac_rho/pi/bac_h);
  clear ind_aob_nob
  %
  % continue with the solid debries consecutive to death of cells
  ind_dead = find(bac_s == 5);
  if ~isempty(ind_dead)
    % radii of the solid debries
    bac_r(ind_dead) = sqrt(bac_e_d(ind_dead)/bac_rhod/pi/bac_h);
  end
  clear ind_dead
  %
  time_gro = toc;
  %
  % Division
  % store the old number of cells (before division)
  tic
  nold = bac_n;
  for i=1:nold
    % excretion of EPS (before division)
    if (bac_s(i)==1 && bac_r(i) > bac_ra(i)*1.4142135)  % excretion of EPS is based on volume/radius, not mass
      % vol EPS = vol active cell
      % choose an angle for the excretion
      fi = rand*pi;
      % increase number of cells with one
      bac_n = bac_n + 1;
      % coordinates (polar), species and radius of the new cell (bac_n)
      bac_x(bac_n)  = bac_x(i) + 2*bac_r(i)*cos(fi);
      bac_y(bac_n)  = bac_y(i) + 2*bac_r(i)*sin(fi);
      bac_s(bac_n)  = 4;
      bac_c(bac_n)  = bac_c(i);
      fraction      = DivisionFraction + 0.2*(rand-0.5);    % fraction of EPS excreted
      bac_m(bac_n)  = 0;                                    % active mass
      bac_e_d(bac_n)  = fraction*bac_e_d(i);                % EPS where is the case (for non producing species is 0)
      bac_ra(bac_n) = 0;                                    % radius of the active biomass
      bac_r(bac_n)  = sqrt(bac_e_d(bac_n)/bac_rhoe/pi/bac_h); % radius of the active biomass + EPS
      %
      % move old parent cell (ng)
      bac_e_d(i) = bac_e_d(i) - bac_e_d(bac_n);
      bac_r(i) = sqrt((bac_m(i)/bac_rho+bac_e_d(i)/bac_rhoe)/pi/bac_h);
    end % EPS excretion
    %
    if (bac_s(i)==1 && bac_m(i) > bac_mmax(1)) || ...   % HET cell division (after mass)
        (bac_s(i)==2 && bac_m(i) > bac_mmax(2)) || ...  % AOB cell division
        (bac_s(i)==3 && bac_m(i) > bac_mmax(3))         % NOB cell division
      % choose an angle for the division
      fi = rand*pi;
      % increase number of cells with one
      bac_n = bac_n + 1;
      % coordinates (polar), species and radius of the new cell (bac_n)
      bac_x(bac_n)  = bac_x(i) + bac_r(i)*cos(fi);
      bac_y(bac_n)  = bac_y(i) + bac_r(i)*sin(fi);
      bac_s(bac_n)  = bac_s(i);
      bac_c(bac_n)  = bac_c(i);
      fraction      = DivisionFraction + 0.2*(rand-0.5);
      bac_m(bac_n)  = fraction*bac_m(i);                   % active mass
      bac_e_d(bac_n)  = fraction*bac_e_d(i);                   % EPS where is the case (for non producing species is 0)
      bac_ra(bac_n) = sqrt(bac_m(bac_n)/pi/bac_rho/bac_h); % radius of the active biomass
      bac_r(bac_n)  = sqrt((bac_m(bac_n)/bac_rho+bac_e_d(bac_n)/bac_rhoe)/pi/bac_h);    % radius of the active biomass + EPS
      % move old parent cell (ng)
      bac_x(i) = bac_x(i) - bac_r(i)*cos(fi); % half distance to the left
      bac_y(i) = bac_y(i) - bac_r(i)*sin(fi); % half distance downward
      %keep the overall mass balanced
      bac_m(i) = bac_m(i) - bac_m(bac_n);
      bac_e_d(i) = bac_e_d(i) - bac_e_d(bac_n);
      % compute the decreased radii (active biomass + EPS)
      bac_ra(i) = sqrt(bac_m(i)/pi/bac_rho/bac_h);
      bac_r(i)  = sqrt((bac_m(i)/bac_rho+bac_e_d(i)/bac_rhoe)/pi/bac_h);
    end % division
    %
  end % all cells were tested
  time_div=toc;
  %
  % calculate center of mass for colonies
  tic
  %
  % clones of debris are marked back with zero
  bac_c(bac_s==5) = 0;
  % calculate center of mass for colonies
  max_c = max(bac_c);
  mtot_c  = zeros(1,max_c); % bac_nc
  mxtot_c = zeros(1,max_c); % bac_nc
  mytot_c = zeros(1,max_c); % bac_nc
  area    = zeros(1,max_c); % bac_nc
  how_many = zeros(1,max_c); % how many cells are in a clone
  for i_c=1:max_c
    ind_c = find(bac_c==i_c);
    n_c = length(ind_c);
    
    bac_x_c = bac_x(ind_c);  bac_y_c = bac_y(ind_c);
    bac_m_c = bac_m(ind_c);  bac_r_c = bac_r(ind_c);
    bac_s_c = bac_s(ind_c);  bac_c_c = bac_c(ind_c);
    
    mtot_c(i_c)  = sum(bac_m_c);
    mxtot_c(i_c) = sum(bac_m_c.*bac_x_c);
    mytot_c(i_c) = sum(bac_m_c.*bac_y_c);
    area(i_c) = pi*sum(bac_r_c.^2);
    how_many(i_c) = n_c;
  end
  %
  xc = mxtot_c./mtot_c;     % 0x coordinate of the center of the virtual suprastucture
  yc = mytot_c./mtot_c;     % 0y coordinate of the center of the virtual suprastucture
  for ihm = 1:length(how_many)
    if how_many(ihm) <= 2
      R_col = Rcol(1);
    else
      R_col = Rcol(2);
    end
    rc(ihm) = sqrt(area(ihm)*R_col/pi);  % radius of the virtual suprastructure; R_col is a factor to include also the voids
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Moving the particulate cells in 2 steps: first inside the colonies of
  % AOB and NOB, and then in the whole floc
  %====================================
  % shoving the cells inside the microcolonies with Java routine
  tic
  iter = 0;
  max_c = max(bac_c);
  for i_c=1:max_c
    ind_c = find(bac_c==i_c);
    n_c = length(ind_c);
    if n_c>1
      bac_x_c = bac_x(ind_c);  bac_y_c = bac_y(ind_c);
      bac_m_c = bac_m(ind_c);  bac_r_c = bac_r(ind_c);
      bac_s_c = bac_s(ind_c);  bac_c_c = bac_c(ind_c);
      Results = QSobj.shoving2Dinside(iterShoving, 1, n_c, bac_x_c, bac_y_c, bac_m_c, bac_r_c, bac_s_c, i_c, xc(i_c), yc(i_c), rc(i_c));
      bac_x(ind_c) = Results.bac_x';    bac_y(ind_c) = Results.bac_y';
      iter  = Results.iter;
    end
    
  end
  
  time_sho1=toc;
  %
  % recalculate center of mass for colonies
  %
  n_collonies = max(bac_c);
  mtot_c  = zeros(1,n_collonies); % bac_nc
  mxtot_c = zeros(1,n_collonies); % bac_nc
  mytot_c = zeros(1,n_collonies); % bac_nc
  area    = zeros(1,n_collonies); % bac_nc
  how_many = zeros(1,n_collonies); % how many cells are in a clone
  %
  ind = [find(bac_s==2) find(bac_s==3)]; % vector of indexes of the aob and nob clones
  while ~isempty(ind) % shrink this vector after working a clone
    jind = find(bac_c(ind(1))==bac_c);
    mtot_c(bac_c(ind(1)))  = sum(bac_m(jind));
    mxtot_c(bac_c(ind(1))) = sum(bac_m(jind).*bac_x(jind));
    mytot_c(bac_c(ind(1))) = sum(bac_m(jind).*bac_y(jind));
    area(bac_c(ind(1))) = pi*sum(bac_r(jind).^2);
    tc(bac_c(ind(1))) = bac_s(ind(1));  % this has the type of the clone
    how_many(bac_c(ind(1))) = length(jind);
    ind(ismember(ind, jind)) = [];      % shrink the vector, to elliminate this collony
  end
  %
  xc = mxtot_c./mtot_c;        % 0x coordinate of the center of the virtual suprastucture
  yc = mytot_c./mtot_c;        % 0y coordinate of the center of the virtual suprastucture
  rc = sqrt(area*R_col/pi);
  for ihm = 1:length(how_many)
    if how_many(ihm) <= 2
      R_col = Rcol(1);
    else
      R_col = Rcol(2);
    end
    rc(ihm) = sqrt(area(ihm)*R_col/pi);  % radius of the virtual suprastructure; 1.5 is a factor to include also the voids
  end
  %
  %   rc = sqrt(area*R_col/pi);
  %
  pos_1 = find(bac_s==1);      % finds the position of the heterotrophs in the vector including all the cells
  pos_4 = find(bac_s==4);      % finds the position of the EPS in the vector including all the cells
  pos_5 = find(bac_s==5);      % finds the position of the solid debries in the vector including all the cells
  n145 = length([pos_1 pos_4 pos_5]); % number of HET + EPS + solid debries individuals
  %
  bac_mshov = bac_m + bac_e_d;   % the mass of particles shoving, considering also eps and debries
  % adds the coordinates for the virtual suprastructure of colonies (AOB and NOB) to the ones of HET and EPS
  bac_xs = [bac_x([pos_1 pos_4 pos_5])  xc];
  % adds the coordinates for the virtual suprastructure of colonies (AOB and NOB) to the ones of HET and EPS
  bac_ys = [bac_y([pos_1 pos_4 pos_5])  yc];
  bac_ms = [bac_mshov([pos_1 pos_4 pos_5]) mtot_c];
  % adds the radius the virtual suprastructure of colonies (AOB and NOB) to the ones of HET and EPS
  bac_rs = [bac_r([pos_1 pos_4 pos_5])  rc];
  bac_t = length(bac_xs);      % length of the vector for transport of cells
  %
  %================================
  % shoving HET and EPS particles with the AOB and NOB microcolonies with Java routine
  %
  tic
  %Results = QSobj.shoving2D_colony(iterShoving, 1, bac_t, bac_xs, bac_ys, bac_ms, bac_rs, 0.2);
  %
  Results = QSobj.shoving2D_colony(iterShoving, 1, bac_t, bac_xs, bac_ys, bac_ms, bac_rs, 1);
  bac_xs = Results.bac_x';    bac_ys = Results.bac_y';    iter_colony  = Results.iter;
  %
  time_sho2=toc;
  %
  bac_x([pos_1 pos_4 pos_5]) = bac_xs(1:n145); % "re-writing" the 0x position in the initial vector, for HET and EPS + solid debries
  bac_y([pos_1 pos_4 pos_5]) = bac_ys(1:n145); % "re-writing" the 0y position in the initial vector, for HET and EPS + solid debries
  %
  bac_xc = bac_xs(n145+1:bac_t) - xc;    % reads from bac_xc the distance with which the clones should move on 0x
  bac_yc = bac_ys(n145+1:bac_t) - yc;    % reads from bac_yc the distance with which the clones should move on 0y
  %
  % moving theAOB and NOB cells in the colonies
  %
  ind = [find(bac_s==2) find(bac_s==3)];
  bac_x(ind) = bac_x(ind) +  bac_xc(bac_c(ind));
  bac_y(ind) = bac_y(ind) +  bac_yc(bac_c(ind));
  %
  clear ind;
  %
  % end of the shoving processes
  %
  %%%%%%%%%%%%%%%%%%%%%%
  % "re-centering" the floc
  %
  %%% center of mass for all the cells
  mtot_floc = sum(bac_m+bac_e_d);
  mxtot_floc = sum((bac_m+bac_e_d).*bac_x);
  mytot_floc = sum((bac_m+bac_e_d).*bac_y);
  %
  xc_flocOld = mxtot_floc/mtot_floc;
  yc_flocOld = mytot_floc/mtot_floc;
  % moving the center of mass for the floc in the center of
  % computational field
  xc_floc = Lx/2;
  yc_floc = Ly/2;
  %
  dxc = xc_flocOld - xc_floc; % the distance on 0x to move all the cells;
  dyc = yc_flocOld - yc_floc; % the distance on 0y to move all the cells;
  
  bac_x = bac_x - dxc;        % these are the new coordinates for the cells of the coming group
  bac_y = bac_y - dyc;
  %
  %
  ind = [find(bac_s==2) find(bac_s==3)];  % vector of indexes of the aob and nob clones
  while ~isempty(ind)                     % shrink this vector after working a clone
    jind = find(bac_c(ind(1))==bac_c);
    mtot_c(bac_c(ind(1)))  = sum(bac_m(jind));
    mxtot_c(bac_c(ind(1))) = sum(bac_m(jind).*bac_x(jind));
    mytot_c(bac_c(ind(1))) = sum(bac_m(jind).*bac_y(jind));
    area(bac_c(ind(1))) = pi*sum(bac_r(jind).^2);
    tc(bac_c(ind(1))) = bac_s(ind(1));  % this has the type of the clone
    how_many(bac_c(ind(1))) = length(jind);
    ind(ismember(ind, jind)) = [];% shrink the vector, to elliminate this collony
  end
  %
  xc = mxtot_c./mtot_c;      % 0x coordinate of the center of the virtual suprastucture
  yc = mytot_c./mtot_c;      % 0y coordinate of the center of the virtual suprastucture
  %
  for ihm = 1:length(how_many)
    if how_many(ihm) <= 2
      R_col = Rcol(1);
    else
      R_col = Rcol(2);
    end
    rc(ihm) = sqrt(area(ihm)*R_col/pi);  % radius of the virtual suprastructure;
  end
  %
  %   rc = sqrt(area*R_col/pi);
  %
  % Put cell masses on a square grid => biomass concentration
  %
  c_bac = zeros(Ny,Nx,5);
  solids = zeros(Ny,Nx);
  for i=1:bac_n                   % loop on all the cells
      jj = floor(bac_x(i)/dx)+1;    % find index of the cell on the grid
      ii = floor(bac_y(i)/dy)+1;
      if bac_s(i) < 4
          c_bac(ii,jj,bac_s(i)) = bac_m(i);
      else
          c_bac(ii,jj,bac_s(i)) = bac_e_d(i); % add cell mass
      end
      %
      solids(ii,jj)=1;
  end
  %
  % Transform biomass in biomass concentration and area in cell fraction
  c_bac = c_bac/(dx*dy*dz);
  %
  solids_2 = solids;
  solids_3 = solids;
  %
  % aici e procesul de dilatare
  for dil=1:10
    for i=2:Nx-1
      for j=2:Ny-1
        if solids_3(i,j)==1
          solids_2(i,j)   = 1;
          solids_2(i+1,j) = 1;
          solids_2(i-1,j) = 1;
          solids_2(i,j+1) = 1;
          solids_2(i,j-1) = 1;
        end
      end
    end
    solids_3 = solids_2;
  end
  %
  % aici punem '2' pentru stratul limita
  for i=1:Nx
    for j=1:Ny
      if solids(i,j)==0 && solids_2(i,j)==1
        solids(i,j) = 2;
      end
    end
  end
  %
  
  for i=1:Nx*Ny  %loop on all the grid elements
    xx = grid(1,i);
    yy = grid(2,i);
    %
    jj = floor(xx/dx)+1;  % find index of the cell on the grid
    ii = floor(yy/dy)+1;
    %
    if solids(ii,jj)== 0
      f_eff(ii,jj) = 1000;
    else
      f_eff(ii,jj) = 1;
    end
  end
  %
  % Recalculate the maximum colony/floc radius; this gives a better besh
  % around the floc, and let coarse mesh outside
  %
  R_floc = max(sqrt((bac_x-Lx/2).^2 + (bac_y-Ly/2).^2)) + 10e-6;  % maximum distance + 10 microns
  %
  for i=1:5
    save_gridded_data(c_bac(:,:,i), xg, yg, PathOutput, ['c_bac_' num2str(i) '_'], iterGrowth+1);
  end
  save_gridded_data(f_eff(:,:), xg, yg, PathOutput, 'f_eff_', iterGrowth+1);
  %
  t = t + dt;
  t_R = [t_R; t];
  %
  % save the data for the reactor concentrations
  file = fopen(nameReactorFile,'a');
  fprintf(file,'%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e',...
    t_R(end), nFlocs(end), c_s_R(end), c_o2_R(end), c_nh4_R(end), c_no2_R(end), c_no3_R(end));
  fprintf(file,'\n');
  fclose(file);
  %
  %
  pHET = [pHET length(find(bac_s == 1)) / length(bac_s)];
  pAOB = [pAOB length(find(bac_s == 2)) / length(bac_s)];
  pNOB = [pNOB length(find(bac_s == 3)) / length(bac_s)];
  pEPS = [pEPS length(find(bac_s == 4)) / length(bac_s)];
  pDEAD = [pDEAD length(find(bac_s == 5)) / length(bac_s)];
  %
  % save the data for the statistics in the floc
  %
  file = fopen(nameFlocStatistics,'a');
  fprintf(file,'%15.5e%15.5e%15.5e%15.5e%15.5e',...
    pHET(end), pAOB(end), pNOB(end), pEPS(end), pDEAD(end));
  fprintf(file,'\n');
  fclose(file);
  %
  %
  % [iterGrowth iter iter_colony bac_n ndetach ndetisland]
  ['iterGrowth=' num2str(iterGrowth) '  innerIter=' num2str(iter) '  outerIter=' num2str(iter_colony)...
    '  totalPartic=' num2str(bac_n) '   detachNo=' num2str(ndetach) '    islandNo=' num2str(ndetisland) ]
  %
  if mod(iterGrowth, ShowImage) == 0
    %
    % Plot solution of cells on the concentration field
    %
    fh=figure('PaperPosition',[0.01 0.0 3.5 3.5],'Position',[ 100 100 690 590 ],'Color','w','PaperPositionMode','auto');
    hax=axes('Position',[ 0.05 0.05 0.85 0.85 ],'Layer','top');
    %
    chem_name = 'c_no2';
    chem_title= 'C_N_O_2';
    %
    min_lim = 0;
    %max_lim = c_no2_reactor_init;
    max_lim = 3*c_no2_reactor_init;
    %
    hp = postplot(fem, ...
      'tridata',{['if(' chem_name '>' num2str(max_lim) ',' num2str(max_lim) ',' chem_name ')'],'cont','internal'}, ...
      'tridlim',[min_lim max_lim], ...
      'trimap','jet(1024)', ...
      'tribar','off', ...
      'solnum','end', ...
      'axis',[0,Lx,0,Ly]);
    %
    axis equal; axis([0 Lx 0 Ly]);
    set(gca,'Clim',[min_lim max_lim]);
    colorbar;
    hold on;
    %
    % indices of the aob and nob bacteria
    %
    bac_nc = max(bac_c);       % number of clones (only for AOB and NOB)
    mtot_c  = zeros(1,bac_nc); %
    mxtot_c = zeros(1,bac_nc); %
    mytot_c = zeros(1,bac_nc); %
    area    = zeros(1,bac_nc); %
    how_many = zeros(1,bac_nc); % how many cells are in a clone
    %
    ind_aob_nob = [find(bac_s==2) find(bac_s==3)];  % vector of indices of the aob and nob clones
    while ~isempty(ind_aob_nob)                     % shrink this vector after working a clone
      jind = find(bac_c(ind_aob_nob(1))== bac_c);   % finds how many cells belong to a cluster; always 1 as ind_aob_nob is shrinking (line 571)
      mtot_c(bac_c(ind_aob_nob(1)))  = sum(bac_m(jind));
      mxtot_c(bac_c(ind_aob_nob(1))) = sum(bac_m(jind).*bac_x(jind));
      mytot_c(bac_c(ind_aob_nob(1))) = sum(bac_m(jind).*bac_y(jind));
      area(bac_c(ind_aob_nob(1)))    = pi*sum(bac_r(jind).^2); %area of the cluster
      tc(bac_c(ind_aob_nob(1))) = bac_s(ind_aob_nob(1)); % this has the type of the clone
      how_many(bac_c(ind_aob_nob(1))) = length(jind);
      ind_aob_nob(ismember(ind_aob_nob, jind)) = [];% shrink the vector, to elliminate this collony
    end
    %
    xc = mxtot_c./mtot_c;      % 0x coordinate of the center of the virtual suprastucture
    yc = mytot_c./mtot_c;      % 0y coordinate of the center of the virtual suprastucture
    for ihm = 1:length(how_many)
      if how_many(ihm) <= 2
        R_col = Rcol(1);
      else
        R_col = Rcol(2);
      end
      rc(ihm) = sqrt(area(ihm)*R_col/pi);  % radius of the virtual suprastructure; R_col is a factor to include also the voids
    end
    %
    for ii=1:bac_nc
      if ~isnan(xc(ii)) && ~isnan(yc(ii))
        if tc(ii)==2
          colcol = [0.8 0.8 1];
        elseif tc(ii)==3
          colcol = [1 0.8 0.8];
        end
        rectangle('Curvature',[1 1],'Position',[xc(ii)-rc(ii) yc(ii)-rc(ii) 2*rc(ii) 2*rc(ii)],'LineWidth',0.1,'EdgeColor','k','FaceColor',colcol);
      end
      hold on;
    end
    %
    for ii=1:bac_n
      rectangle('Curvature',[1 1],'Position',[bac_x(ii)-bac_r(ii) bac_y(ii)-bac_r(ii) 2*bac_r(ii) 2*bac_r(ii)],'LineWidth',1,'EdgeColor','none','FaceColor',col{4});
      
      if bac_s(ii)<4
        rectangle('Curvature',[1 1],'Position',[bac_x(ii)-bac_ra(ii) bac_y(ii)-bac_ra(ii) 2*bac_ra(ii) 2*bac_ra(ii)],...
          'LineWidth',0.1,'EdgeColor','none','FaceColor',col{bac_s(ii)});    %*(1 - bac_c(i)*0.01));
      end
      % inert representation
      if bac_s(ii)==5
        rectangle('Curvature',[1 1],'Position',[bac_x(ii)-bac_r(ii) bac_y(ii)-bac_r(ii) 2*bac_r(ii) 2*bac_r(ii)],...
          'LineWidth',1,'EdgeColor','none','FaceColor',col{5});
      end
      hold on;
    end
    %
    axis equal; axis([0 Lx 0 Ly]);
    set(gca,'XTick',[0:50e-6:Lx],'FontSize',14);
    set(gca,'YTick',[0:50e-6:Ly],'FontSize',14);
    %   set(gca,'XTickLabel',{'0';'50';'100';'150';'200'});
    %   set(gca,'YTickLabel',{'0';'50';'100';'150';'200'});
    set(gca,'visible','on','Layer','top','Box','on');
    %
    ht = title(['Bacteria at t = ' num2str(t,'%5.3f') ' days']);
    set(ht,'FontSize',16);
    hx = xlabel('x [\mum]');  set(hx,'FontSize',16);
    hy = ylabel('y [\mum]');  set(hy,'FontSize',16);
    %
    drawnow;
    %
    % saves the plot of cells in the concentration field at every second iteration
    %
    print(gcf, '-dpng', '-r600', [PathOutput '\cells.' num2str(iterGrowth,'%04d') '.png']);
    %       print(gcf, '-depsc', [PathOutput '\cells.' num2str(iterGrowth,'%04d') '.eps']);
    %
    hold off;
    close all
    %
    iplot=1;
    plot_1;
  end
  %
  %
end
% here ends the time integration
% ================================================================================
figure
plot(t_R, c_s_R, '-r')
hold on
plot(t_R, c_o2_R, '-b')
plot(t_R, c_nh4_R, ':g')
plot(t_R, c_no2_R, '-.g')
plot(t_R, c_no3_R, '-.g')
xlabel('time, day')
ylabel('Reactor concentration, kg/m^3')
legend('S', 'O_2', 'NH_4', 'NO_2', 'NO_3')
%
figure
semilogy(t_R, nFlocs, '-r')
%
t
