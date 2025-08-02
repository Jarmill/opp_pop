%test the system from https://tobiasgeyer.org/RaKG23_3LOPPsforLCLfilters.pdf

%this is a given LCL filter
%not a custom filter

%grid frequency 
f0 = 50;
w0 = f0*2*pi;

%filter
Lf = 0.35e-3;
Rf = 0.3e-3;
C = 420e-6;
Rc = 4e-3;
Lt = 526.41e-6;
Rt = 16.54e-3;

%grid/machine
Lg = 349.19e-6;
Rg = 10.97e-3;

%compound parameters
Lgt = Lg + Lt;
R1 = Rf + Rc;
R2 = Rc + Rg + Rt;

%% design the plant

%grid-connected converter
Arg = [-R1, Rc, -1, 0, 0;
    Rc, -R2, 1, -1, 0;
    1, -1, 0, 0, 0;
    0, 0, 0, 0, -w0;
    0, 0, 0, w0, 0;];
Asg = diag(1./[Lf, Lgt, C, 1, 1]);

Ag = Arg*Asg;
Bg = [1; 0; 0; 0; 0];
Cg = [0, 1, 0, 0, 0];
Dg = 0;

sysg = ss(Ag, Bg, Cg, Dg);


%no grid connection
Ar = [-R1, Rc, -1;
    Rc, -R2, 1;
    1, -1, 0];
As = diag(1./[Lf, Lgt, C]);

A = Ar*As;
B = [1; 0; 0];
C = [0, 1, 0];
D = 0;

sys = ss(A, B, C, D);

%grid-connected converter (in input form
Bgi = [1 0; 0 1; 0 0];
Cg = [0, 1, 0];
Dgi = [0, 0];

sysgi = ss(A, Bgi, Cg, Dgi);

%only damping
Ad = [-Rg/Lg];
Bd = 1;
Cd = 1;
Dd = 0;
sysd = ss(Ad, Bd, Cd, Dd);


%% apply a signal
NP = 6;
NT = 1000;
t = linspace(0, NP/f0, NT)';
u = sin(w0*t);
[yf, to, xo] = lsim(sys, u, t);
[ydf, to, xo] = lsim(sysd, u, t);


%% plots
figure(1)
clf
hold on
bode(sys)
bode(sysg)
bode(sysd)

figure(2)
clf
hold on
plot(t, u);
plot(t, yf);
plot(t, ydf);