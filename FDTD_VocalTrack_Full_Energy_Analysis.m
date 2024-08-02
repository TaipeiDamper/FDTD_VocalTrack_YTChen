%   Start Date 2024/08
%
%   File: FDTD_VocalTrack_Full_Energy_Analysis
%   Auther : Yung-Tao Chen
%
%   the code is developed through the following order
%   (1) just a cylinder
%   (2) Area variation (dA/dx ~= 0)
%   (3) radiation condition at lip term (right end)
%   (4) cross-section area oscillation (not used at right end grid since it is radiation pt)
%   (5) add nasal part
%   (6) collect radiation boundary energy -> failed, fix radiation approx
%   later
%   (7) collect wall loss
%
%
%   %%%%%%% comment %%%%%%%
%   This is the 1-Dimensional Finite Difference Time Domain(FDTD)
%   vocal track simulation.
%   The vocal tract is treated as multiple tubes connected to each other,
%   with wall oscillation(loss), radiation loss and nasal coupling.
%   All energy stored in the system is calculated, and energy conservation
%   is recovered.
%
%   In short, the energy in the system can be calculated in to the
%   following parts.
%   1. E_p: the energy stored in pressure term
%   2. E_v: the energy stored in velocity term
%   3. E_rad_loss: the energy loss at boundary in time (should be power),
%                  controlled by "a1"
%   4. E_rad_osc: the energy stored at the boundary (as an inductance).
%                  controlled by "a2"
%   5. E_m_w: the energy stored at wall mass, controlled by m_w
%   6. E_b_w: the energy loss at wall damper in time (should be power), 
%                 controlled by b_w
%   7. E_k_w: the energy stored at wall stiffness, controlled by k_w
%
%   To set the system lossless, b_w and a1 should be 0 (b_w ~ line 100
%                                                        a1 ~ line 220)
%   To see the inpulse-response, set the input to be delta function (~line 170)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------
%         grid spacing
%        |←    X    →| 
%        |           |   q-0.5   q   q+0.5  q+1   
%        |           |  
%        |           |      update value
%        |           |           ↓
%        p           p           p           p             n+1    -------
%                                ⇑                                   ↑ 
%                           ⇗        ⇖                                 
%    v         v           v           v            v      n+0.5     T
%                                ↑                                        
%                                ↑                                   ↓
%        p           p           p           p              n     -------
%                                                                   time
%
%    v         v           v           v            v      n-0.5  
%
%   glottis                                     boundary

% initialize
clear
close all
clc

% parameters
SR = 44100;             % sampleRate
T = 1/SR;               % time step (sec)
Tf = 2;                 % simulation time [sec]
Nf = floor(SR*Tf);      % number of time steps
t = (0:Nf-1)'/SR;       % (1:SR*Tf)'/SR; time vector
rho = 1.293;            % density of air kg/m^3


c = 340;                % wave speed (sound speed in air) m/s
L = 0.175;              % vocal track length (m)
r = c/L;
f = 100;                % fundamental frequency of the input
Zc = rho*c;
Zc_n = rho*c;
%(5)
L_n = 0.105;            % nasal tract length (m)
nasal_start = 0.08;     % the connecting point position (m), 0 = at glottis
n_ratio = 0.5;          % nasal ratio, if = 0, all air go to lip
                        %              if = 1, all air go to nose



% wall param(4)         % this term for vocal tract wall oscillation
                        % b would be the loss term
                        % if b = 0, the whole scheme should be lossless
m_w = 1.5;              % wall mass/unit area
b_w = 1400;             % wall damping/unit area
% b_w = 0;                % uncomment to set the wall to be lossless
k_w = 10^-5;            % wall stiffness/unit area          
                        % (k_w value only effect the system when b = 0 and -k is really large)

% -----------------------------------------
%%%%%%%%%%%%%%%%%%%%% calculate grids in the system %%%%%%%%%%%%%%%%%%%%%%%%%
% courant number: Lambda = c*k/h = c*T/X and Lamdba -> 1(-) for small
% numerical dispersion(close to 1) and prevent explosive growth(<1)
% c*T/X <= 1  -> c*T <= X 
X = c*T;                % exact grid spacing
N = floor(L/X);         % number or grid
Xnew = L/N;             % the grid spacing actually used in the scheme

% (5)
N_n = floor(L_n/Xnew);                      % number of grid in nasal cavity
connect_grid = floor(nasal_start/L*N +1);   % calculate the connecting grid
                                            % this grid would be the p grid
                                            %                     <- v_n(1)
                                            % v -> p(connect_grid)<- v(connect_grid+1)


% stability check
if (c*T/Xnew > 1)
disp('courant number > 1')                  % show error message if >1
end
disp(c*T/Xnew)

% initialize vectors
p = zeros(1,N);     % pressure grid used in the scheme
v = zeros(1,N);     % particle velocity used in the scheme
y = zeros(Nf,1);    % output vector
% (5)
p_n = zeros(1,N_n); % pressure grid used in nasal cavity
v_n = zeros(1,N_n); % particle velocity used in nasal cavity
y_n = zeros(Nf,1);
y_total = zeros(Nf,1);

% energy analysis
E_v = zeros(Nf,1);      % energy of velocity term
E_p = zeros(Nf,1);      % energy of pressure term
E_total = zeros(Nf,1);  % total energy
% (6)
E_total_timeDifference = zeros(Nf,1);
E_system = zeros(Nf,1);
% (5)
E_v_n = zeros(Nf,1);    % energy of velocity term in nasal cavity 
E_p_n = zeros(Nf,1);    % energy of pressure term in nasal cavity


% (6)
E_rad_loss = zeros(Nf,1);
E_rad_loss_last = 0;
E_rad_osc = zeros(Nf,1);
% (5)
E_rad_loss_n = zeros(Nf,1);
E_rad_osc_n = zeros(Nf,1);
E_rad_loss_last_n = 0;

% (7)
E_wall_loss = zeros(Nf,1);
E_wall_loss_n = zeros(Nf,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the input energy would be sent to the glottis end 

% input at glottis
% f = f.*2.^(0.0012*sin(2*pi*5*t));   % create vibration in input signal, 
% f = f.*2.^((t/0.25)/12);


v_Gl = sin(2*pi*t.*f);
v_Gl = 0.5*(abs(v_Gl) +v_Gl) ;      % this would be the semi-sinwave
%%%%%%% uncomment this to set see input to be impulse %%%%%%%%
% v_Gl = [1;zeros(Nf-1,1)];    % turn on this for impulse input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Area term (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vocal tract profile, non-dimensional [pos S] pairs
% the vocal tract area is calculated in this part

S0 = 0.00025;                               % vocal tract surface area, left end (m^2)

% area data from NSS
% pick vowel needed
% /E/
S = [0 1;0.09 0.4;0.11 2.4;0.24 2.4;0.26 3.2;0.29 3.2;0.32 4.2;...
0.41 4.2;0.47 3.2;0.59 1.8;0.65 1.6;0.71 1.6;0.74 1;0.76 0.8;...
0.82 0.8;0.88 2;0.91 2;0.94 3.2;1 3.2];

% /A/
% S = [0 1;0.03 0.60;0.09 0.4;0.12 1.6;0.18 0.6;0.29 0.2;0.35 0.4;...
% 0.41 0.8;0.47 1;0.50 0.6;0.59 2;0.65 3.2;0.85 3.2;0.94 2;1 2];


% area grid for p(pressure) and v(group velocity)
grid = 0:1/(2*N-1):1;                    % total number of grid used in vocal tract
A = S0*interp1(S(:,1),S(:,2),grid);      % vocal tract area data


% Area term for nasal cavity (5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the nasal area is calculated in this part

S_n = [0 0.5;0.5 0.5;1 0.6;2 1.5;2.7 1.7;3.5 1.7;4 1.5;4.4 2;6 1.7;6.8 2;7.5 1.7;8 2.2;8.5 2;9 2;9.3 1.6;9.7 1.8;10 1.5;10.5 0.2];
% 0 is nose end, need to be inverted
S_n(:,1) = (S_n(end,1)-S_n(:,1)); % invert the order, so the nose end would be 1(in next step)
S_n(:,2) = S_n(:,2)/10000; % convert to cm^2
S_n = flip(S_n,1);


% area grid for p(pressure) and v(group velocity)
grid_n = 0:1/(2*N_n-1):1;                   % total number of grid used in nasal tract
A_n = interp1(S_n(:,1),S_n(:,2),grid_n);    % nasal tract area data


% radiation boundary condition param (3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Aend_r = sqrt(A(end)/pi);                   % take the radias at the lips end of the tube
% (5)
Aend_r_n = sqrt(A_n(end)/pi);               % take the radias at the nose end of the tube

% radiaiton condition parameter
% (roughly) following the equation: v_end+1 = a1*p_end + a2*m_end, dm/dt = p
a1 = 1/(4*0.6133^2);
a2 = c/(0.6133*Aend_r);

a1_n = 1/(4*0.6133^2);
a2_n = c/(0.6133*Aend_r_n);


%%%%%% uncomment to set the radiation end to be lossless %%%%%
% a1 = 0;
% a1_n = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "alpha" is calculated to simplify the equation in the dsp loop
alpha = a1/Zc+a2/Zc*T/2;
% (5)
alpha_n = a1_n/Zc_n+a2_n/Zc_n*T/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dsp loop initialize ------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update p first, then update v in a loop %
% p = p(nt,(n+0.5)k)                      %
% v = v((n+0.5)t, nk)                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% previous loop param (1)
v_last = zeros(1,length(v));        % store data in last loop
p_last = zeros(1,length(p));
% radiation param (3)
m = 0;       % dm/dt = p                          
m_last = 0;
% (5) copy anything above
v_last_n = zeros(1,length(v_n));    % store data in last loop
p_last_n = zeros(1,length(p_n));
% radiation param (3)
m_n = 0;     % dm/dt = p                           
m_last_n = 0;


% Area function (cross-section) (2)
Area_v = A(1:2:end-1);                                              % get area at v grid
Area_v_j_v = Area_v(connect_grid+1)*(1-n_ratio);                    % recalculate junction Area, only for calculate area_p_n
Area_p = (Area_v(2:end) + Area_v(1:end-1))/2;                       % get area from average area at v grid (for energy conservation)
Area_p(connect_grid+1) = (Area_v(connect_grid+2) + Area_v_j_v)/2;   % the area ratio should only change the grid on the right

Area_p = horzcat(Area_p,Area_v(end));                               % since there is no v(n+1)grid, use v(n) area to be the area at p(n)
Area_v_plus = horzcat(Area_v,Area_p(end));                          % add the last value for easier code reading

% (5) copy anything above
% Area function (cross-section) (2)
Area_v_n = A_n(1:2:end-1);                                          % get area at v grid
Area_v_n(1) = Area_v(connect_grid)*n_ratio;                         % recalculate junction Area, only for calculate area_p_n
Area_p_n = (Area_v_n(2:end) + Area_v_n(1:end-1))/2;                 % get area from average area at v grid (for energy conservation)

Area_p_n = horzcat(Area_p_n,Area_v_n(end));                         % since there is no v(n+1)grid, use v(n) area to be the area at p(n)
Area_v_plus_n = horzcat(Area_v_n,Area_p_n(end));                    % add the last value for easier code reading

% Area oscillation term (4)
peri_p = 2*pi*sqrt(Area_p(1:end-1)/pi);                             % perimeter of area at p grid
disp_0 = zeros(1, length(Area_p(1:end-1)));                         % displacement now
disp_f = disp_0;                                                    % displacement in the future
disp_1 = disp_0;                                                    % displacement in the past
disp_2 = disp_0;                                                    % further past
% (5) copy anything above
% Area oscillation term (4)
peri_p_n = 2*pi*sqrt(Area_p_n(1:end-1)/pi);                         % perimeter of area at p grid
disp_0_n = zeros(1, length(Area_p_n(1:end-1)));                     % displacement now
disp_f_n = disp_0_n;                                                % displacement in the future
disp_1_n = disp_0_n;                                                % displacement in the past
disp_2_n = disp_0_n;                                                % further past

% (5) junction that connect all tubes (most important part in nasal coupling)
v_j = 0;  % this term would take the area at area_v(connect_grid+1)
          % then the volume velocity here would be seperate to
          % v(connect_grid+1) and v_n(1), while
          % this is a virtual grid, so this grid would be excluded from
          % energy analysis
v_j_last = 0;


%%%%%%%%%%%%%%%%%%%%%%%%% dsp main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
for n = 1:Nf


% update p -----------------------------------------------------------------
var_a_term = (disp_0-disp_1)/T.*peri_p;  % dA/dt
var_a_term_n = (disp_0_n-disp_1_n)/T.*peri_p_n;  % dA/dt

p(1:end-1) = p_last(1:end-1) - T*rho*c^2/Xnew*(v(2:end).*Area_v(2:end) - v(1:end-1).*Area_v(1:end-1))./Area_p(1:end-1) + (-1)*var_a_term*rho*c^2*T./Area_p(1:end-1);
%                                                                                                                           |-----------   area oscillator term (4) --------|
% (5) deal with the connecting point for the vocal tract and the nasal tract 
% juicy part!
p(connect_grid) = p_last(connect_grid) - T*rho*c^2/Xnew*(v_j.*Area_v(connect_grid+1)- v(connect_grid).*Area_v(connect_grid))./Area_p(connect_grid) + (-1)*var_a_term(connect_grid)*rho*c^2*T./Area_p(connect_grid);
p(connect_grid+1) = p_last(connect_grid+1) - T*rho*c^2/Xnew*(v(connect_grid+2).*Area_v(connect_grid+2)- v_j*(1-n_ratio).*Area_v(connect_grid+1))./Area_p(connect_grid+1) + (-1)*var_a_term(connect_grid+1)*rho*c^2*T./Area_p(connect_grid+1);

% (5) same as p above for nasal tarct
p_n(1:end-1) = p_last_n(1:end-1) - T*rho*c^2/Xnew*(v_n(2:end).*Area_v_n(2:end) - v_n(1:end-1).*Area_v_n(1:end-1))./Area_p_n(1:end-1) + (-1)*var_a_term_n*rho*c^2*T./Area_p_n(1:end-1);
% particularly deal with p_n(1)
p_n(1) = p_last_n(1) - T*rho*c^2/Xnew*(v_n(2).*Area_v_n(2) - v_j*n_ratio.*Area_v(connect_grid+1))./Area_p_n(1) + (-1)*var_a_term_n(1)*rho*c^2*T./Area_p_n(1);


%%%    %%% below part for radiation condition (3) %%%   %%%
% radiation p end
%////% p(end) = p_last(end)-T*rho*c^2/X*(1*( p_last(end)*(alpha)+a2/Zc*(2*m_last) ) - v(end)*(1+1)); // same as below (rearrange)
p(end) = p_last(end)*(1-0.5*T*rho*c^2/Xnew*alpha)-T*rho*c^2/Xnew/Area_p(end)*(Area_p(end)*a2/Zc*(m_last)- v(end)*(Area_v(end)));
p(end) = p(end)/(1+0.5*T*rho*c^2/Xnew/1*alpha);

p_n(end) = p_last_n(end)*(1-0.5*T*rho*c^2/Xnew*alpha_n)-T*rho*c^2/Xnew/Area_p_n(end)*(Area_p_n(end)*a2_n/Zc_n*(m_last_n)- v_n(end)*(Area_v_n(end)));
p_n(end) = p_n(end)/(1+0.5*T*rho*c^2/Xnew/1*alpha_n);

%%%%% % lossless lip end (closed) (replace/uncomment the value above if needed)%%%%%%%
% % p(end) = p_last(end)- T*rho*c^2/Xnew*(0 - v(end)*Area_v(end))/Area_p(end) + (-1)*var_a_term(end)*rho*c^2*T/Area_p(end);
% p(end) = p_last(end)- T*rho*c^2/Xnew*(0 - v(end)*Area_v(end))/Area_p(end);
% 
% % (5)
% % p_n(end) = p_last_n(end)- T*rho*c^2/Xnew*(0 - v_n(end)*Area_v_n(end))/Area_p_n(end) + (-1)*var_a_term_n(end)*rho*c^2*T/Area_p_n(end);
% p_n(end) = p_last_n(end)- T*rho*c^2/Xnew*(0 - v_n(end)*Area_v_n(end))/Area_p_n(end);

% lossless lip end (opened)
% p(end) = 0;
% p_n(end) = 0;

%%%%%  %%% above part for radiation condition %%%   %%%%%%%%%%%%%%%%%%%%%%%%

% update m (3)
m = m_last + T/2*(p(end) + p_last(end));
% (5)
m_n = m_last_n + T/2*(p_n(end) + p_last_n(end));



% update variation area (4)
disp_f = ((2*m_w/T^2+b_w/T-k_w)*disp_0 - m_w/T^2*disp_1 + p(1:end-1))/(m_w/T^2+b_w/T);
% (5)
disp_f_n = ((2*m_w/T^2+b_w/T-k_w)*disp_0_n - m_w/T^2*disp_1_n + p_n(1:end-1))/(m_w/T^2+b_w/T);

% energy for pressure term (1)
E_p(n) = 1/2/rho/c^2*sum((p(1:connect_grid).^2).*(Area_v_plus(2:connect_grid+1)+Area_v_plus(1:connect_grid))/2)*Xnew; 
E_p(n) = E_p(n) + 1/2/rho/c^2*sum((p(connect_grid+2:end).^2).*(Area_v_plus(connect_grid+3:end)+Area_v_plus(connect_grid+2:end-1))/2)*Xnew; 
E_p(n) = E_p(n) + 1/2/rho/c^2*sum((p(connect_grid+1).^2).*(Area_v_plus(connect_grid+2)+Area_v_j_v)/2)*Xnew; 
% (6)
E_rad_osc(n) = a2/Zc/2*Area_p(end)*m^2;                 % energy stored at the boundary grid, no "X" term(not in sigma/boundary)

if (n>1)
E_rad_loss(n-1) =  a1/Zc*Area_p(end)*((p(end)+p_last(end))/2)^2;
end

% (5)
E_p_n(n) = 1/2/rho/c^2*sum((p_n.^2).*(Area_v_plus_n(2:end)+Area_v_plus_n(1:end-1))/2)*Xnew; 
% (6)
E_rad_osc_n(n) = a2_n/Zc_n/2*Area_p_n(end)*m_n^2;       % energy stored at the boundary grid, no "X" term(not in sigma/boundary)

if (n>1)
E_rad_loss_n(n-1) = a1_n/Zc_n*Area_p_n(end)*((p_n(end)+p_last_n(end))/2)^2;
end

% energy for wall oscillator term (4)
E_wall_m = sum( peri_p*m_w/(2*T^2).*(disp_f.*disp_0 + disp_1.*disp_0 - disp_0.*disp_0 - disp_f.*disp_1) )*Xnew; % energy for m
E_wall_k = sum( peri_p*k_w/2.*(disp_0.^2) )*Xnew;       % energy for k
% (5)
E_wall_m_n = sum( peri_p_n*m_w/(2*T^2).*(disp_f_n.*disp_0_n + disp_1_n.*disp_0_n - disp_0_n.*disp_0_n - disp_f_n.*disp_1_n) )*Xnew; % energy for m
E_wall_k_n = sum( peri_p_n*k_w/2.*(disp_0_n.^2) )*Xnew; % energy for k

%(7)
if(n>1)
E_wall_loss(n-1) = sum(peri_p*b_w/2/T^2.*((disp_f-disp_1).*(disp_0-disp_1)))*Xnew;
E_wall_loss_n(n-1) = sum(peri_p_n*b_w/2/T^2.*((disp_f_n-disp_1_n).*(disp_0_n-disp_1_n)))*Xnew;
end

% update v ----------------------------------------------------------------- 
v(2:end) = v_last(2:end) - T/rho/Xnew*(p(2:end)-p(1:end-1)); 
% (5)
v_n(2:end) = v_last_n(2:end) - T/rho/Xnew*(p_n(2:end)-p_n(1:end-1)); 

% another juicy part!
v(connect_grid+1) = v_last(connect_grid+1) - T/rho/Xnew*(p(connect_grid+1)-p(connect_grid)); % next velocity grid from the junction in vocal tract
% (5)
v_n(1) = v_last_n(1) - T/rho/Xnew*(p_n(1)-p(connect_grid));     % first grid in nasal cavity
% update v_j from the combination of v(connect_grid+1) and v_n(1) (5)
v_j = v(connect_grid+1)*(1-n_ratio) + v_n(1)*n_ratio;

% input at glottis
v(1) = v_Gl(n);     % take input v (when input = inpulse = lossless closed end)

% calculate energy for velocity term
E_v(n) = rho/2*sum(v(1:connect_grid).*v_last(1:connect_grid).*Area_v(1:connect_grid))*Xnew + rho/2*sum(v(connect_grid+2:end).*v_last(connect_grid+2:end).*Area_v(connect_grid+2:end))*Xnew ;  % energy for velocity grid
E_v(n) = E_v(n) + rho/2*(v_j*v_j_last*Area_v(connect_grid+1))*Xnew; % add energy at connect grid
% (5)
E_v_n(n) = rho/2*sum(v_n(2:end).*v_last_n(2:end).*Area_v_n(2:end))*Xnew;  % energy for velocity grid




%%%%%%%%%% only keep this to remind myself that it existed before, not used now %%%%%%%%%%
% v(end) = v_last(end)-T/rho/X*(0-p(end));   % open end (this would not be used since the right end would be p)
% v(end) = v_last(end)+a1/rho*(p(end)-p_last(end))+T*a2/rho*p_last(end);  % radiation end
%%%%%%%%%% only keep this to remind myself that it existed before, not used now %%%%%%%%%%

% calculate energy total in this step
E_total(n) = E_p(n) + E_v(n)+ E_wall_k + E_wall_m;
% (5)
E_total(n) = E_total(n) + E_p_n(n) + E_v_n(n)+ E_wall_k_n + E_wall_m_n;
E_total(n) = E_total(n) + E_rad_osc(n) + E_rad_osc_n(n);
% (6)
if (n>1)
E_total_timeDifference(n-1) = (E_total(n)-E_total(n-1))/T;
E_system(n-1) = E_total_timeDifference(n-1) + E_rad_loss(n-1) + E_rad_loss_n(n-1);
E_system(n-1) = E_system(n-1) + E_wall_loss(n-1) + E_wall_loss_n(n-1);
end
% pass variable for next loop --------------------------------------------
% (1)
v_last = v;         
p_last = p;
m_last = m;

% (4)
disp_2 = disp_1;
disp_1 = disp_0;
disp_0 = disp_f;

% (5) copy
% (1)
v_last_n = v_n;         
p_last_n = p_n;
m_last_n = m_n;

% (4)
disp_2_n = disp_1_n;
disp_1_n = disp_0_n;
disp_0_n = disp_f_n;

% (5)
v_j_last = v_j;

% (6)
E_rad_loss_last = E_rad_loss(n);
E_rad_loss_last_n = E_rad_loss_n(n);

% take the output at lip pressure grid
y(n) = p(end);
y_n(n) = p_n(end);
y_total(n) = p(end) + p_n(end);
% y_total(n) =  p_n(end);
end
% loop end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
figure(1)
subplot(3,1,1)
plot(t,v_Gl,'k')
xlim([0 0.1])
xlabel('time (sec)')
ylabel('Amplitude')
title('input signal')

subplot(3,1,2);
plot(t,y,'r')
xlim([0 0.1])
xlabel('time (sec)')
ylabel('Amplitude')
title('output signal at lips')

subplot(3,1,3);
plot(t,y_n,'r')
xlim([0 0.1])
xlabel('time (sec)')
ylabel('Amplitude')
title('output signal at nose')

figure(2)
subplot(3,1,1);
in_freq = abs(fft(v_Gl));
out_freq = abs(fft(y));
out_freq_n = abs(fft(y_n));
out_freq_all = abs(fft(y+y_n));
f_vector = (0:Nf-1)*SR/Nf;
plot(f_vector,10*log10(in_freq),'k')
xlabel('frequency (Hz)')
ylabel('Amplitude')
title('input signal')

subplot(3,1,2);
plot(f_vector,10*log10(out_freq),'r')
% semilogx(f_vector,10*log10(out_freq),'r')
xlabel('frequency (Hz)')
ylabel('Amplitude')
title('output signal at lips')
xlim([0 5000])

subplot(3,1,3);
% plot(f_vector,10*log10(out_freq_all),'r')
plot(f_vector,10*log10(out_freq_n),'r')
% semilogx(f_vector,10*log10(out_freq_n),'r')
xlabel('frequency (Hz)')
ylabel('Amplitude')
title('output signal at nose')
xlim([0 5000])

figure(3)
subplot(2,1,1)
plot(t(2:end),E_total(2:end),'k')
xlabel('time')
ylabel('total Energy in body')
title('output signal')

subplot(2,1,2)
plot(t(2:end),(E_total(2:end)-E_total(2))/E_total(2))
xlim([0 0.1])
xlabel('time')
ylabel('variation of total Energy in body')
title('output signal')

figure(4)
subplot(2,1,1)
plot(t(2:end-1),E_system(2:end-1),'k')
xlabel('time')
ylabel('total Energy')
title('total Energy in system (including loss)')

subplot(2,1,2)
plot(t(3:end-1),(E_system(3:end-1)-E_system(2))/E_system(2))
xlim([0 0.1])
xlabel('time')
ylabel('variation of total Energy in system')
title('total Energy variation in system (including loss)')

soundsc(y_total,SR);
% soundsc(y,SR)
n_ratio