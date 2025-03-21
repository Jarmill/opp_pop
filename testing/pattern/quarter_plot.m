   % alpha[1] = 0.6055937525415025
   % alpha[2] = 0.6213040492559058
   % alpha[3] = 0.6370187386620668
   % alpha[4] = 0.6527302458124219
   % alpha[5] = 0.6684416046082732
   % alpha[6] = 0.6841567523061838
   % alpha[7] = 0.6998670059248276
   % alpha[8] = 1.5629421252042202
   %    du[1] = 1.0
   %    du[2] = 0.0
   %    du[3] = 1.0
   %    du[4] = 0.0
   %    du[5] = 1.0
   %    du[6] = 0.0
   %    du[7] = 1.0
   %    du[8] = 0.0

   % a = 0.6674572160283855;
   % alpha = [a, pi-a, a+pi, 2*pi-a];
   % u = [0 1 0 -1 0];

   % a = [0.6055937525415025;
   %     0.6213040492559058
   % 0.6370187386620668;
   % 0.6527302458124219;
   % 0.6684416046082732;
   % 0.6841567523061838;
   % 0.6998670059248276;
   % 1.5629421252042202];
   % u0 = 0;
   % uv = [0 1 0 1 0 1 0 1 0]';

   %modulation = 1;
      % a = [0.6514286491352599; 0.6671438104191006; 0.6828588859056817];
% uv = [0 1 0 1]';
% u0 = 0;

modulation = 0.55;
 a = [0.5764188652596557;
     1.4940129440240595;
     1.5098427393429221;
     1.5256472497294145;
     1.562892065118792]
      uv = [0 1 0 1 0 1]';

   aa_quad = [0; kron(a, [1; 1]); pi/2];
    u_quad = kron(uv, [1; 1]);
    aa_all = [aa_quad; pi - (aa_quad(end:-1:1)); pi + aa_quad; 2*pi - aa_quad(end:-1:1)];
    u_all = [u_quad; u_quad(end:-1:1); -u_quad; -u_quad(end:-1:1)];

   % [na, nb] = pulse_harmonics(3, u, alpha);

   figure(34)
   clf
   hold on
   plot(aa_all, u_all, 'linewidth', 1)
   xlim([0, 2*pi])

   theta = linspace(0, 2*pi, 500);
   plot(theta, modulation*sin(theta), ':k', 'linewidth', 2);
   % xlabel('\theta')
   % ylabel('u(\theta)')
   axis off

   figure(40)
   clf
   
   hold on
   plot3(cos(aa_all), sin(aa_all), u_all)
   plot3(cos(theta), sin(theta), modulation*sin(theta), ':k', 'linewidth', 2)
   axis off
   xlim([-1, 1])
   ylim([-1, 1])
   % zlim([-1, 1])
   % pbaspect([1 1 1])
   view(3)
   
