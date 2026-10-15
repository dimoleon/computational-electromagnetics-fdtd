c0 = 299792458;
eps0 = 8.85418781762e-12;
mu0 = 1.256637061436e-06;

freq = 3.75e9;
cutoff_freq = 3e9;

lambda0 = c0 / freq;
b = c0 / (2 * cutoff_freq);
a = 6 * lambda0;

dx = lambda0 / 10;
nx = ceil(a / dx);
dx = a / nx;
ny = ceil(b / dx);
dy = b / ny;

x_int = linspace(0.0, a, nx + 1);
y_int = linspace(0.0, b, ny + 1);

dt_courant = 1.0 / (c0 * sqrt(1/dx^2 + 1/dy^2));

source_x_idx = nx / 2 + 1;                      % f matlab
source_profile = sin(pi * y_int / b);
