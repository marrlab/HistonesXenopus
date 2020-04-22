function [model] = histonesXenopusmock_MM_0_d1d2d3_r1r2r3_syms()

% create state syms
syms 'K20P' 'K20m1' 'K20m2' 'K20m3'
% create state vector
model.sym.x = [K20P K20m1 K20m2 K20m3];

% create parameter syms
syms 'a' 'b' 'd' 'K20m1_0' 'K20m2_0' 'K20m3_0' 'noise' 'd1_m' 'd2_m' 'd3_m' 'r1_m' 'r2_m' 'r3_m'  
% create parameter vector
model.sym.p = [a b d K20m1_0 K20m2_0 K20m3_0 noise d1_m d2_m d3_m r1_m r2_m r3_m];
% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10';

% create symbolic variable for time
syms 't'
model.sym.xdot = sym(zeros(size(model.sym.x)));
% piecewise defined function
model.sym.xdot(1) = -r1_m*K20P+d1_m*K20m1+log(2)/(a+(b*t)/(d+t))...
    *(K20P+K20m1+K20m2+K20m3)-(log(2)*(a*(d+t)^2+b*t^2)/(a*(d+t)+b*t)^2)*K20P;
model.sym.xdot(2) = r1_m*K20P-(r2_m+d1_m)*K20m1+d2_m*K20m2...
    -(log(2)*(a*(d+t)^2+b*t^2)/(a*(d+t)+b*t)^2)*K20m1;
model.sym.xdot(3) = r2_m*K20m1-(r3_m+d2_m)*K20m2+d3_m*K20m3...
    -(log(2)*(a*(d+t)^2+b*t^2)/(a*(d+t)+b*t)^2)*K20m2;
model.sym.xdot(4) = r3_m*K20m2-d3_m*K20m3...
    -(log(2)*(a*(d+t)^2+b*t^2)/(a*(d+t)+b*t)^2)*K20m3;
 
%create initial conditions
model.sym.x0 = sym(zeros(size(model.sym.x)));
model.sym.x0(1) = 0.1/(0.1+K20m1_0+K20m2_0+K20m3_0);
model.sym.x0(2) = K20m1_0/(0.1+K20m1_0+K20m2_0+K20m3_0);
model.sym.x0(3) = K20m2_0/(0.1+K20m1_0+K20m2_0+K20m3_0);
model.sym.x0(4) = K20m3_0/(0.1+K20m1_0+K20m2_0+K20m3_0);

model.sym.k = [];

%observables (relative)
model.sym.y = sym(zeros(1,1));
model.sym.y(1) = log(K20P)-log(K20P+K20m1+K20m2+K20m3);
model.sym.y(2) = log(K20m1)-log(K20P+K20m1+K20m2+K20m3);
model.sym.y(3) = log(K20m2)-log(K20P+K20m1+K20m2+K20m3);
model.sym.y(4) = log(K20m3)-log(K20P+K20m1+K20m2+K20m3);


model.sym.sigma_y = noise*sym(ones(4,1));


end
