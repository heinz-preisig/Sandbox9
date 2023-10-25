function output_vars = main
  % Index sets
  N = cell(1, 6);
  N(1)  = 1;
  N(2)  = 2;
  N(3)  = 3;
  N(4)  = 4;
  N(5)  = 5;
  N(6)  = 6;
  N_lbl = "N";
  N_x_S = cell(1, 6);
  N_x_S(1)  = [1];
  N_x_S(2)  = [1];
  N_x_S(3)  = [1];
  N_x_S(4)  = [1];
  N_x_S(5)  = [1];
  N_x_S(6)  = [1];
  N_x_S_lbl = "N_x_S";
  A = cell(1, 5);
  A(1)  = 1;
  A(2)  = 2;
  A(3)  = 3;
  A(4)  = 4;
  A(5)  = 5;
  A_lbl = "A";
  A_x_S = cell(1, 5);
  A_x_S(1)  = [1];
  A_x_S(2)  = [1];
  A_x_S(3)  = [1];
  A_x_S(4)  = [1];
  A_x_S(5)  = [1];
  A_x_S_lbl = "A_x_S";
  S = cell(1, 1);
  S(1)  = 1;
  S_lbl = "S";

  N_cons = [1, 2];
  N_dyna = [3, 4, 5, 6];
  A_diff = [1, 2, 3, 4, 5];

  % Variables
  % starting time
  to = 0.0;

  % initial species
  no = MultiDimVar({N_x_S_lbl}, {N_x_S});
  no(6) = [0.0];
  no(5) = [0.0];
  no(3) = [0.0];
  no(4) = [0.0];

  % end time
  te = 100.0;

  % link variable nProd to interface reactions >>> macroscopic
  _nProd = MultiDimVar({N_x_S_lbl}, {N_x_S});
  _nProd(4) = [0.0];
  _nProd(5) = [0.0];
  _nProd(3) = [0.0];
  _nProd(6) = [0.0];

  % species related incidence matrix
  F_NS_AS = MultiDimVar({N_x_S_lbl, A_x_S_lbl}, {N_x_S, A_x_S});
  F_NS_AS.value(3, 1) = -1;
  F_NS_AS.value(4, 1) = 1;
  F_NS_AS.value(4, 2) = -1;
  F_NS_AS.value(5, 2) = 1;
  F_NS_AS.value(5, 3) = -1;
  F_NS_AS.value(6, 3) = 1;
  F_NS_AS.value(1, 4) = -1;
  F_NS_AS.value(3, 4) = 1;
  F_NS_AS.value(6, 5) = -1;
  F_NS_AS.value(2, 5) = 1;

  % link variable k_d_Fick to interface material >>> macroscopic
  _k_d_Fick = MultiDimVar({N_x_S_lbl}, {N_x_S});
  _k_d_Fick(5) = [6e-06];
  _k_d_Fick(4) = [6e-06];
  _k_d_Fick(2) = [6e-06];
  _k_d_Fick(1) = [6e-06];
  _k_d_Fick(3) = [6e-06];

  % difference operator for species topology
  D_NS_AS = MultiDimVar({N_x_S_lbl, A_x_S_lbl}, {N_x_S, A_x_S});
  D_NS_AS.value(3, 1) = -1;
  D_NS_AS.value(4, 1) = 1;
  D_NS_AS.value(4, 2) = -1;
  D_NS_AS.value(5, 2) = 1;
  D_NS_AS.value(5, 3) = -1;
  D_NS_AS.value(6, 3) = 1;
  D_NS_AS.value(1, 4) = -1;
  D_NS_AS.value(3, 4) = 1;
  D_NS_AS.value(6, 5) = -1;
  D_NS_AS.value(2, 5) = 1;

  % molar composition
  c = MultiDimVar({N_x_S_lbl}, {N_x_S});
  c(2) = [0.0];
  c(1) = [0.01];

  % cross sectional area yz
  Ayz = MultiDimVar({N_lbl}, {N});
  Ayz(5) = [1e-07];
  Ayz(4) = [1e-07];
  Ayz(2) = [1e-07];
  Ayz(1) = [1e-07];
  Ayz(3) = [1e-07];

  % fundamental state -- volume
  V = MultiDimVar({N_lbl}, {N});
  V(4) = [1e-06];
  V(5) = [1e-06];
  V(3) = [1e-06];
  V(6) = [1e-06];

  % fundamental state -- molar mass
  n = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % differential species balance
  dndt = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % net molar convectional mass flow
  fnc = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % net diffusional mass flow
  fnd = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % diffusional mass flow in a given stream
  fnd_AS = MultiDimVar({A_x_S_lbl}, {A_x_S});


% Integrators
  % Initial conditions
  phi_0(1:4) = (no(N_dyna)).value;

  % Integration interval
  integration_interval = [to, te];
  dt = 0.1;
  options = odeset('InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-3, 'AbsTol', 1e-3);
  % Integrator routine
  [t, phi] = ode45(@f, integration_interval, phi_0, options);

  % Integrand
  function dphidt = f(t, phi) 
    n(N_dyna) = phi(1:4);

    fnc(N_dyna) = E_75(F_NS_AS, fnc_AS, N_dyna, A_diff);

    c(N_dyna) = E_44(n, V, N_dyna);

    fnd_AS(A_diff) = E_152(_k_d_Fick, D_NS_AS, c, Ayz, indexunion(N_dyna, N_cons), A_diff);

    fnd(N_dyna) = E_69(F_NS_AS, fnd_AS, N_dyna, A_diff);

    dndt(N_dyna) = E_76(fnc, fnd, _nProd, N_dyna);

    fprintf("%g    ", t);
    fprintf("%g    ", c.value(N_dyna));
    fprintf("\n");
    dphidt(1:4) = (dndt(N_dyna)).value;
  endfunction
endfunction

% Functions for the equations
function sol = E_75(F_NS_AS, fnc_AS, N, A)
  sol = reduceproduct(F_NS_AS(N, A), "A_x_S", fnc_AS(A));
endfunction

function sol = E_44(n, V, N)
  sol = khatrirao(1 ./ V(N), n(N));
endfunction

function sol = E_152(_k_d_Fick, D_NS_AS, c, Ayz, N, A)
  sol = reduceproduct(khatrirao(Ayz(N), -1 .* _k_d_Fick(N)) .* D_NS_AS(N, A), "N_x_S", c(N));
endfunction

function sol = E_69(F_NS_AS, fnd_AS, N, A)
  sol = reduceproduct(F_NS_AS(N, A), "A_x_S", fnd_AS(A));
endfunction

function sol = E_76(fnc, fnd, _nProd, N)
  sol = fnc(N) + fnd(N) + _nProd(N);
endfunction


% Auxiliar functions
function result = indexunion(varargin)
  result = varargin{1};
  for i = 2:nargin
    result = union(result, varargin{i});
  endfor
endfunction