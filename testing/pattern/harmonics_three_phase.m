N = 9000;
th = linspace(0, 2*pi, N);
thh = th(1:N/2);
thq = th(1:N/4);
CM = true;

% % x = sin(th) + 0.5*cos(3*th);

% f = @(t) sin(t) - 0.5*sin(5*t) - 0.2*sin(7*t);
% f = @(t) sin(t) - 0.5*sin(5*t) - 0.2*sin(11*t);
if ~CM
    f = @(t) 2*sin(t) - sin(5*t) - 0.2*sin(11*t);
else
    f = @(t) sin(t) - 0.5*sin(5*t) - 0.2*sin(11*t) + 0.2*sin(9*t);
end

x = f(th);

xa = x;
xb = f(th - 2*pi/3);
xc = f(th - 4*pi/3);
% xb = circshift(x, N/3);
% xc = circshift(x, 2*N/3);


% index = 0:5;

n = 1;
sa = sin(n*th);
ca = cos(n*th);
sb = sin(n*(th+2*pi/3));
cb = cos(n*(th+2*pi/3));
sc = sin(n*(th-2*pi/3));
cc = cos(n*(th-2*pi/3));

%% full-wave index shifts
a_base = simps(th, sa .* xa)/pi;
a_base_b = simps(th, sa .* xb)/pi;
a_b = simps(th, sb .* xa)/pi;
a_base_c = simps(th, sa .* xc)/pi;
a_c = simps(th, sc .* xa)/pi;
b_base_b = simps(th, ca .* xb)/pi;
b_b = simps(th, cb .* xa)/pi;
b_c = simps(th, cc .* xa)/pi;


%% half-wave index shifts

a_base_h = 2*simps(th(1:N/2), sa(1:N/2) .* xa(1:N/2))/pi;
% a_base_b_h = 2*simps(th, sa(1:N/2) .* xb(1:N/2))/pi;
a_bh = 2*simps(th(1:N/2), sb(1:N/2) .* xa(1:N/2))/pi;
a_ch = 2*simps(th(1:N/2), sc(1:N/2) .* xa(1:N/2))/pi;
b_bh = 2*simps(th(1:N/2), cb(1:N/2) .* xa(1:N/2))/pi;
b_ch = 2*simps(th(1:N/2), cc(1:N/2) .* xa(1:N/2))/pi;


%% quarter-wave index shifts

a_bq = 4*simps(th(1:N/4), sb(1:N/4) .* xa(1:N/4))/pi;
% a_bh = simps(th, sa .* xa)/pi;
% b_base_b = simps(th, ca .* xb)/pi;
% b_b = simps(th, cb .* xa)/pi;
