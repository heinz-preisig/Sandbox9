function output_vars = main
  % Index sets
  N = cell(1, 2);
  N(1)  = 1;
  N(2)  = 2;
  N_lbl = "N";
  N_x_S = cell(1, 2);
  N_x_S(1)  = [1];
  N_x_S(2)  = [1];
  N_x_S_lbl = "N_x_S";
  A = cell(1, 1);
  A(1)  = 1;
  A_lbl = "A";
  A_x_S = cell(1, 1);
  A_x_S(1)  = [1];
  A_x_S_lbl = "A_x_S";
  S = cell(1, 1);
  S(1)  = 1;
  S_lbl = "S";

  N_cons = [1];
  N_dyna = [2];
  A_diff = [1];

  % Variables
  % starting time
  to = 0.0;

  % end time
  te = 100.0;

  % initial species
  no = MultiDimVar({N_x_S_lbl}, {N_x_S});
  no(2) = [0.0];

  % link variable nProd to interface reactions >>> macroscopic
  _nProd = MultiDimVar({N_x_S_lbl}, {N_x_S});
  _nProd(2) = [0.0];

  % net diffusional mass flow
  fnd = MultiDimVar({N_x_S_lbl}, {N_x_S});
  fnd(2) = [0.0];

  % net molar convectional mass flow
  fnc = MultiDimVar({N_x_S_lbl}, {N_x_S});
  fnc(2) = [2.0];

  % fundamental state -- molar mass
  n = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % differential species balance
  dndt = MultiDimVar({N_x_S_lbl}, {N_x_S});


% Integrators
  % Initial conditions
  phi_0(1:1) = (no(N_dyna)).value;

  % Integration interval
  integration_interval = [to, te];
  dt = 0.1;
  options = odeset('InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-3, 'AbsTol', 1e-3);
  % Integrator routine
  [t, phi] = ode45(@f, integration_interval, phi_0, options);

  % Integrand
  function dphidt = f(t, phi) 
    % simple 
     dndt(N_dyna) = E_76(_nProd, fnd, fnc, N_dyna);

    fprintf("%g    ", t);
    fprintf("%g    ", dndt.value(N_dyna));
    fprintf("%g    ", phi);
    fprintf("\n");
    dphidt(1:1) = (dndt(N_dyna)).value;
  endfunction
endfunction

% Functions for the equations
    % simple 
function sol = E_76(_nProd, fnd, fnc, N)
  sol = plus(plus(fnc, fnd), _nProd);
endfunction


% Auxiliar functions
function result = indexunion(varargin)
  result = varargin{1};
  for i = 2:nargin
    result = union(result, varargin{i});
  endfor
endfunction