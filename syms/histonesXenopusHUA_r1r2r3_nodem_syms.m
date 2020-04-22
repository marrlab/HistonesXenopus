function [model] = histonesXenopusHUA_r1r2r3_nodem_syms()

% create state syms
syms 'K20P' 'K20m1' 'K20m2' 'K20m3'
% create state vector
model.sym.x = [K20P K20m1 K20m2 K20m3];

% create parameter syms
syms 'noise' 'r1_h' 'r2_h' 'r3_h' 
% create parameter vector
model.sym.p = [noise r1_h r2_h r3_h];
% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10';

% create symbolic variable for time
syms 't'
model.sym.xdot = sym(zeros(size(model.sym.x)));
% piecewise defined function
model.sym.xdot(1) = -r1_h*K20P;
model.sym.xdot(2) = r1_h*K20P-(r2_h)*K20m1;
model.sym.xdot(3) = r2_h*K20m1-(r3_h)*K20m2;
model.sym.xdot(4) = r3_h*K20m2;

syms 'k1' 'k2' 'k3' 'k4'
model.sym.k = [k1, k2, k3, k4];

model.sym.x0 = [k1, k2, k3, k4];

%observables (relative)
model.sym.y = sym(zeros(1,1));
model.sym.y(1) = log(K20P)-log(K20P+K20m1+K20m2+K20m3);
model.sym.y(2) = log(K20m1)-log(K20P+K20m1+K20m2+K20m3);
model.sym.y(3) = log(K20m2)-log(K20P+K20m1+K20m2+K20m3);
model.sym.y(4) = log(K20m3)-log(K20P+K20m1+K20m2+K20m3);


model.sym.sigma_y = noise*sym(ones(4,1));


end
