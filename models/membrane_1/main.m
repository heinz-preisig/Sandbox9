function output_vars = main
  % Index sets
  N = cell(1, 10);
  N(1)  = 1;
  N(2)  = 2;
  N(3)  = 3;
  N(4)  = 4;
  N(5)  = 5;
  N(6)  = 6;
  N(7)  = 7;
  N(8)  = 8;
  N(9)  = 9;
  N(10)  = 10;
  N_lbl = "N";
  N_x_S = cell(1, 10);
  N_x_S(1)  = [1];
  N_x_S(2)  = [1];
  N_x_S(3)  = [1];
  N_x_S(4)  = [1];
  N_x_S(5)  = [1];
  N_x_S(6)  = [1];
  N_x_S(7)  = [1];
  N_x_S(8)  = [1];
  N_x_S(9)  = [1];
  N_x_S(10)  = [1];
  N_x_S_lbl = "N_x_S";
  A = cell(1, 9);
  A(1)  = 1;
  A(2)  = 2;
  A(3)  = 3;
  A(4)  = 4;
  A(5)  = 5;
  A(6)  = 6;
  A(7)  = 7;
  A(8)  = 8;
  A(9)  = 9;
  A_lbl = "A";
  A_x_S = cell(1, 9);
  A_x_S(1)  = [1];
  A_x_S(2)  = [1];
  A_x_S(3)  = [1];
  A_x_S(4)  = [1];
  A_x_S(5)  = [1];
  A_x_S(6)  = [1];
  A_x_S(7)  = [1];
  A_x_S(8)  = [1];
  A_x_S(9)  = [1];
  A_x_S_lbl = "A_x_S";
  S = cell(1, 1);
  S(1)  = 1;
  S_lbl = "S";

  N_dyna = [1, 2, 7, 8, 9, 10];
  N_cons = [3, 4, 5, 6];
  A_diff = [1, 2, 3, 4, 5];
  A_conv = [6, 7, 8, 9];

  % Variables
  % initial species
  no = MultiDimVar({N_x_S_lbl}, {N_x_S});
  no(9) = [0.0];
  no(2) = [0.0];
  no(8) = [0.0];
  no(7) = [0.0];
  no(10) = [0.0];
  no(1) = [0.0];

  % starting time
  to = 0.0;

  % end time
  te = 100.0;

  % link variable nProd to interface reactions >>> macroscopic
  _nProd = MultiDimVar({N_x_S_lbl}, {N_x_S});
  _nProd(9) = [0.0];
  _nProd(2) = [0.0];
  _nProd(8) = [0.0];
  _nProd(10) = [0.0];
  _nProd(7) = [0.0];
  _nProd(1) = [0.0];

  % species related incidence matrix
  F_NS_AS = MultiDimVar({N_x_S_lbl, A_x_S_lbl}, {N_x_S, A_x_S});
  F_NS_AS.value(7, 1) = -1;
  F_NS_AS.value(8, 1) = 1;
  F_NS_AS.value(8, 2) = -1;
  F_NS_AS.value(9, 2) = 1;
  F_NS_AS.value(9, 3) = -1;
  F_NS_AS.value(10, 3) = 1;
  F_NS_AS.value(1, 4) = -1;
  F_NS_AS.value(7, 4) = 1;
  F_NS_AS.value(10, 5) = -1;
  F_NS_AS.value(2, 5) = 1;

  % cross sectional area yz
  Ayz = MultiDimVar({N_lbl}, {N});
  Ayz(5) = [1e-07];
  Ayz(4) = [1e-07];
  Ayz(1) = [1e-07];
  Ayz(3) = [1e-07];
  Ayz(2) = [1e-07];
  Ayz(7) = [1e-07];
  Ayz(9) = [1e-07];
  Ayz(8) = [1e-07];
  Ayz(6) = [1e-07];

  % difference operator for species topology
  D_NS_AS = MultiDimVar({N_x_S_lbl, A_x_S_lbl}, {N_x_S, A_x_S});
  D_NS_AS.value(7, 1) = -1;
  D_NS_AS.value(8, 1) = 1;
  D_NS_AS.value(8, 2) = -1;
  D_NS_AS.value(9, 2) = 1;
  D_NS_AS.value(9, 3) = -1;
  D_NS_AS.value(10, 3) = 1;
  D_NS_AS.value(1, 4) = -1;
  D_NS_AS.value(7, 4) = 1;
  D_NS_AS.value(10, 5) = -1;
  D_NS_AS.value(2, 5) = 1;

  % link variable k_d_Fick to interface material >>> macroscopic
  _k_d_Fick = MultiDimVar({N_x_S_lbl}, {N_x_S});
  _k_d_Fick(5) = [6e-06];
  _k_d_Fick(4) = [6e-06];
  _k_d_Fick(1) = [6e-06];
  _k_d_Fick(3) = [6e-06];
  _k_d_Fick(2) = [6e-06];

  % fundamental state -- volume
  V = MultiDimVar({N_lbl}, {N});
  V(9) = [1e-06];
  V(2) = [0.001];
  V(8) = [1e-06];
  V(7) = [1e-06];
  V(10) = [1e-06];
  V(1) = [0.001];

  % numerical value one half
  onehalf = [0.5];

  % molar composition
  c = MultiDimVar({N_x_S_lbl}, {N_x_S});
  c(4) = [0.0];
  c(6) = [0.0];
  c(5) = [0.0];
  c(3) = [0.01];

  % link variable density to interface material >>> macroscopic
  _density = MultiDimVar({N_lbl}, {N});
  _density(7) = [1000.0];
  _density(9) = [1000.0];
  _density(8) = [1000.0];
  _density(6) = [1000.0];

  % link variable kc_x to interface material >>> macroscopic
  _kc_x = MultiDimVar({N_lbl}, {N});
  _kc_x(7) = [0.001];
  _kc_x(9) = [0.001];
  _kc_x(8) = [0.001];
  _kc_x(6) = [0.001];

  % fundamental state -- molar mass
  n = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % differential species balance
  dndt = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % net diffusional mass flow
  fnd = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % net molar convectional mass flow
  fnc = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % diffusional mass flow in a given stream
  fnd_AS = MultiDimVar({A_x_S_lbl}, {A_x_S});

  % molar convectional mass flow in the given stream
  fnc_AS = MultiDimVar({A_x_S_lbl}, {A_x_S});

  % concentration in convectional flow
  c_AS = MultiDimVar({A_x_S_lbl}, {A_x_S});

  % volumetric flow
  fV = MultiDimVar({A_lbl}, {A});

  % flow direction of convectional flow
  d = MultiDimVar({A_lbl}, {A});


% Integrators
  % Initial conditions
  phi_0(1:6) = (no(N_dyna)).value;

  % Integration interval
  integration_interval = [to, te];
  dt = 0.1;
  options = odeset('InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-3, 'AbsTol', 1e-3);
  % Integrator routine
  [t, phi] = ode45(@f, integration_interval, phi_0, options);

  % Integrand
  function dphidt = f(t, phi) 
    n(N_dyna) = phi(1:6);

    c(N_dyna) = E_44(n, V, N_dyna);

    fV(A_conv) = E_67(Ayz, _density, p, _kc_x, D, indexunion(N_cons, N_dyna), A_conv);

    d(A_conv) = E_72(F_N_A, p, indexunion(N_cons, N_dyna), A_conv);

    fnd_AS(A_diff) = E_152(Ayz, D_NS_AS, c, _k_d_Fick, N_dyna, A_diff);

    c_AS(A_conv) = E_73(onehalf, c, d, F_NS_AS, indexunion(N_cons, N_dyna), A_conv);

    fnd(N_dyna) = E_69(fnd_AS, F_NS_AS, N_dyna, indexunion(A_conv, A_diff));

    fnc_AS(A_conv) = E_74(c_AS, fV, A_conv);

    fnc(N_dyna) = E_75(fnc_AS, F_NS_AS, N_dyna, indexunion(A_conv, A_diff));

    dndt(N_dyna) = E_76(_nProd, fnd, fnc, N_dyna);

    fprintf("%g    ", t);
    fprintf("%g    ", c.value(N_dyna));
    fprintf("\n");
    dphidt(1:6) = (dndt(N_dyna)).value;
  endfunction
endfunction

% Functions for the equations
function sol = E_44(n, V, N)
  sol = khatrirao(1 ./ V(N), n(N));
endfunction

function sol = E_67(Ayz, _density, p, _kc_x, D, N, A)
  sol = reduceproduct(1 ./ _density(N) .* _kc_x(N) .* Ayz(N) .* D(N, A), "N", p(N));
endfunction

function sol = E_72(F_N_A, p, N, A)
  sol = sign(reduceproduct(F_N_A(N, A), "N", p(N)));
endfunction

function sol = E_152(Ayz, D_NS_AS, c, _k_d_Fick, N, A)
  sol = reduceproduct(khatrirao(Ayz(N), -1 .* _k_d_Fick(N)) .* D_NS_AS(N, A), "N_x_S", c(N));
endfunction

function sol = E_73(onehalf, c, d, F_NS_AS, N, A)
  sol = reduceproduct((onehalf .* (F_NS_AS(N, A) - khatrirao(d(A), abs(F_NS_AS(N, A))))), "N_x_S", c(N));
endfunction

function sol = E_69(fnd_AS, F_NS_AS, N, A)
  sol = reduceproduct(F_NS_AS(N, A), "A_x_S", fnd_AS(A));
endfunction

function sol = E_74(c_AS, fV, A)
  sol = khatrirao(fV(A), c_AS(A));
endfunction

function sol = E_75(fnc_AS, F_NS_AS, N, A)
  sol = reduceproduct(F_NS_AS(N, A), "A_x_S", fnc_AS(A));
endfunction

function sol = E_76(_nProd, fnd, fnc, N)
  sol = fnc(N) + fnd(N) + _nProd(N);
endfunction


% Auxiliar functions
function result = indexunion(varargin)
  result = varargin{1};
  for i = 2:nargin
    result = union(result, varargin{i});
  endfor
endfunction