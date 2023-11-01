function output_vars = main
  % Index sets
  N = cell(1, 12);
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
  N(11)  = 11;
  N(12)  = 12;
  N_lbl = "N";
  N_x_S = cell(1, 12);
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
  N_x_S(11)  = [1];
  N_x_S(12)  = [1];
  N_x_S_lbl = "N_x_S";
  A = cell(1, 13);
  A(1)  = 1;
  A(2)  = 2;
  A(3)  = 3;
  A(4)  = 4;
  A(5)  = 5;
  A(6)  = 6;
  A(7)  = 7;
  A(8)  = 8;
  A(9)  = 9;
  A(10)  = 10;
  A(11)  = 11;
  A(12)  = 12;
  A(13)  = 13;
  A_lbl = "A";
  A_x_S = cell(1, 13);
  A_x_S(1)  = [1];
  A_x_S(2)  = [1];
  A_x_S(3)  = [1];
  A_x_S(4)  = [1];
  A_x_S(5)  = [1];
  A_x_S(6)  = [1];
  A_x_S(7)  = [1];
  A_x_S(8)  = [1];
  A_x_S(9)  = [1];
  A_x_S(10)  = [1];
  A_x_S(11)  = [1];
  A_x_S(12)  = [1];
  A_x_S(13)  = [1];
  A_x_S_lbl = "A_x_S";
  S = cell(1, 1);
  S(1)  = 1;
  S_lbl = "S";

  N_cons = [1, 2];
  N_dyna = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
  A_diff = [1, 2, 3, 5, 6, 7, 8, 11, 12, 13];
  A_conv = [4, 9, 10];

  % Variables
  % initial species
  no = MultiDimVar({N_x_S_lbl}, {N_x_S});
  no(8) = [0.0];
  no(11) = [0.0];
  no(6) = [1.0];
  no(7) = [0.0];
  no(12) = [1.0];
  no(4) = [0.0];
  no(5) = [0.0];
  no(3) = [0.0];
  no(10) = [1.0];
  no(9) = [0.0];

  % starting time
  to = 0.0;

  % end time
  te = 20.0;

  % link variable nProd to interface reactions >>> macroscopic
  _nProd = MultiDimVar({N_x_S_lbl}, {N_x_S});
  _nProd(8) = [0.0];
  _nProd(11) = [0.0];
  _nProd(6) = [0.0];
  _nProd(7) = [0.0];
  _nProd(12) = [0.0];
  _nProd(4) = [0.0];
  _nProd(5) = [0.0];
  _nProd(3) = [0.0];
  _nProd(10) = [0.0];
  _nProd(9) = [0.0];

  % net molar convectional mass flow
  fnc = MultiDimVar({N_x_S_lbl}, {N_x_S});
  fnc(8) = [0.0];
  fnc(11) = [0.0];
  fnc(6) = [0.0];
  fnc(7) = [0.0];
  fnc(12) = [0.0];
  fnc(4) = [0.0];
  fnc(5) = [0.0];
  fnc(3) = [0.0];
  fnc(10) = [0.0];
  fnc(9) = [0.0];

  % species related incidence matrix
  F_NS_AS = MultiDimVar({N_x_S_lbl, A_x_S_lbl}, {N_x_S, A_x_S});
  F_NS_AS.value(3, 1) = -1;
  F_NS_AS.value(4, 1) = 1;
  F_NS_AS.value(4, 2) = -1;
  F_NS_AS.value(5, 2) = 1;
  F_NS_AS.value(6, 3) = -1;
  F_NS_AS.value(3, 3) = 1;
  F_NS_AS.value(7, 5) = -1;
  F_NS_AS.value(8, 5) = 1;
  F_NS_AS.value(8, 6) = -1;
  F_NS_AS.value(9, 6) = 1;
  F_NS_AS.value(10, 7) = -1;
  F_NS_AS.value(7, 7) = 1;
  F_NS_AS.value(12, 8) = -1;
  F_NS_AS.value(11, 8) = 1;
  F_NS_AS.value(5, 11) = -1;
  F_NS_AS.value(8, 11) = 1;
  F_NS_AS.value(9, 12) = -1;
  F_NS_AS.value(11, 12) = 1;
  F_NS_AS.value(11, 13) = -1;
  F_NS_AS.value(2, 13) = 1;

  % difference operator for species topology
  D_NS_AS = MultiDimVar({N_x_S_lbl, A_x_S_lbl}, {N_x_S, A_x_S});
  D_NS_AS.value(3, 1) = -1;
  D_NS_AS.value(4, 1) = 1;
  D_NS_AS.value(4, 2) = -1;
  D_NS_AS.value(5, 2) = 1;
  D_NS_AS.value(6, 3) = -1;
  D_NS_AS.value(3, 3) = 1;
  D_NS_AS.value(7, 5) = -1;
  D_NS_AS.value(8, 5) = 1;
  D_NS_AS.value(8, 6) = -1;
  D_NS_AS.value(9, 6) = 1;
  D_NS_AS.value(10, 7) = -1;
  D_NS_AS.value(7, 7) = 1;
  D_NS_AS.value(12, 8) = -1;
  D_NS_AS.value(11, 8) = 1;
  D_NS_AS.value(5, 11) = -1;
  D_NS_AS.value(8, 11) = 1;
  D_NS_AS.value(9, 12) = -1;
  D_NS_AS.value(11, 12) = 1;
  D_NS_AS.value(11, 13) = -1;
  D_NS_AS.value(2, 13) = 1;

  % link variable k_d_Fick_A to interface material >>> macroscopic
  _k_d_Fick_A = MultiDimVar({A_x_S_lbl}, {A_x_S});
  _k_d_Fick_A(12) = [1.0];
  _k_d_Fick_A(2) = [1.0];
  _k_d_Fick_A(3) = [1.0];
  _k_d_Fick_A(13) = [1.0];
  _k_d_Fick_A(1) = [1.0];
  _k_d_Fick_A(5) = [1.0];
  _k_d_Fick_A(6) = [1.0];
  _k_d_Fick_A(11) = [1.0];
  _k_d_Fick_A(8) = [1.0];
  _k_d_Fick_A(7) = [1.0];

  % molar composition
  c = MultiDimVar({N_x_S_lbl}, {N_x_S});
  c(2) = [0.0];

  % cross sectional area yz (arc)
  Ayz_A = MultiDimVar({A_lbl}, {A});
  Ayz_A(12) = [1.0];
  Ayz_A(2) = [1.0];
  Ayz_A(3) = [1.0];
  Ayz_A(13) = [5.0];
  Ayz_A(1) = [1.0];
  Ayz_A(5) = [1.0];
  Ayz_A(6) = [1.0];
  Ayz_A(11) = [1.0];
  Ayz_A(8) = [1.0];
  Ayz_A(7) = [1.0];

  % fundamental state -- volume
  V = MultiDimVar({N_lbl}, {N});
  V(11) = [1.0];
  V(8) = [1.0];
  V(6) = [1.0];
  V(7) = [1.0];
  V(12) = [1.0];
  V(4) = [1.0];
  V(5) = [1.0];
  V(3) = [1.0];
  V(10) = [1.0];
  V(9) = [1.0];

  % fundamental state -- molar mass
  n = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % differential species balance
  dndt = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % net diffusional mass flow
  fnd = MultiDimVar({N_x_S_lbl}, {N_x_S});

  % diffusional mass flow in a given stream
  fnd_AS = MultiDimVar({A_x_S_lbl}, {A_x_S});


% Integrators
  % Initial conditions
  phi_0(1:10) = (no(N_dyna)).value;

  % Integration interval
  integration_interval = [to, te];
  dt = 1;
  options = odeset('InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-3, 'AbsTol', 1e-3);
  % Integrator routine
  [t, phi] = ode45(@f, integration_interval, phi_0, options);
  save out.txt t phi

  % Integrand
  function dphidt = f(t, phi) 
    n(N_dyna) = phi(1:10);

    c(N_dyna) = E_44(V, n, N_dyna);

    fnd_AS(A_diff) = E_152(D_NS_AS, _k_d_Fick_A, c, Ayz_A, indexunion(N_cons, N_dyna), A_diff);

    fnd(N_dyna) = E_69(fnd_AS, F_NS_AS, N_dyna, indexunion(A_diff, A_conv));

    dndt(N_dyna) = E_76(_nProd, fnd, fnc, N_dyna);

    %fprintf("%g    ", t);
    %fprintf("%g    ", c.value(N_dyna));
    %fprintf("\n");
    dphidt(1:10) = (dndt(N_dyna)).value;
  endfunction
endfunction

% Functions for the equations
function sol = E_44(V, n, N)
  sol = khatrirao(1 ./ V(N), n(N));
endfunction

function sol = E_152(D_NS_AS, _k_d_Fick_A, c, Ayz_A, N, A)
  sol = reduceproduct(khatrirao(Ayz_A(A), -1 .* _k_d_Fick_A(A)) .* D_NS_AS(N, A), "N_x_S", c(N));
endfunction

function sol = E_69(fnd_AS, F_NS_AS, N, A)
  sol = reduceproduct(F_NS_AS(N, A), "A_x_S", fnd_AS(A));
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