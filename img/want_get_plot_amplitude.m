
%output from load_3_level
u =cell(4, 1);
a = cell(4, 1);
M = (0.25:0.25:1)';

u{1} = [0	1	0	1	0	1	0	1	0	-1	0	-1	0	-1	0	-1	0]';
a{1} = [0.740097233420496	0.846562109632576	1.26040344404746	1.38480300858407	1.75678964500572	1.88118920954233	2.29503054395722	2.40149542016930	3.88168988701029	3.98815476322237	4.40199609763725	4.52639566217386	4.89838229859552	5.02278186313213	5.43662319754701	5.54308807375909]';
u{2} = [0	1	0	1	0	1	0	1	0	-1	0	-1	0	-1	0	-1	0]';
a{2} = [0.649098039545527	0.850589287928731	1.16892972309406	1.43484838648996	1.70674426709983	1.97266293049573	2.29100336566106	2.49249461404427	3.79069069313532	3.99218194151852	4.31052237668386	4.57644104007976	4.84833692068962	5.11425558408552	5.43259601925086	5.63408726763406]';
u{3} = [0	1	0	1	0	1	0	1	0	-1	0	-1	0	-1	0	-1	0]';
a{3} = [0.547516875842327	0.808024504451902	1.03641005352050	1.48754112732841	1.65405152626139	2.10518260006929	2.33356814913789	2.59407577774747	3.68910952943212	3.94961715804170	4.17800270711030	4.62913378091820	4.79564417985118	5.24677525365908	5.47516080272768	5.73566843133726]';
u{4} = [0	1	0	1	0	1	0	1	0	-1	0	-1	0	-1	0	-1	0]';
a{4} = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	1.59076253334333	2.33497693624565	2.49530338574537	2.71826955294643	3.56491575423316	3.78788192143422	3.94820837093393	4.69242277383625	4.73235518693313	5.47656958983545	5.63689603933516	5.85986220653622]';

modulation = 1;
th = linspace(0, 2*pi, 1200);
xdes = M*sin(th);
ides = -M*cos(th);

%% compute the current
outq = cell(4, 1);
for i = 1:4
    outq{i} = pulse_current_voltage(u{i}, a{i}, 0);
    outq{i}.uu = kron(u{i}, [1; 1]);    
    outq{i}.aa = [0; kron(a{i}, [1; 1]); 2*pi];

end


%% plot only the quarter-wave
figure(2)
clf
tiledlayout(2, 2)

c = linspecer(4);
% c = c([1, 2,4, 3], :);

c = zeros(4, 3);
c(2, :) = [0.4416    0.7490    0.4322];
c(4, :) = [192, 29, 208]/255;

nexttile
hold on
for i =2:2:4
    plot(th, xdes(i, :), 'k',  'LineWidth', 2, 'Color', c(i, :))
end
xlim([0, 2*pi])
ylabel('$u^*(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);



% c = c(5:end, :);
nexttile
hold on
% plot([pi, pi], [-1, 1], ':k','LineWidth',2)
% plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
% plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
for i = 4:-2:2
    plot(outq{i}.aa, outq{i}.uu, 'LineWidth',3, 'Color', c(i, :))
end
xlim([0, 2*pi])
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);

% axis off
ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);

nexttile
hold on
for i =2:2:4
    plot(th, ides(i, :), 'k',  'LineWidth', 2, 'Color', c(i, :))
end
xlim([0, 2*pi])
ylabel('$I^*(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);

nexttile
hold on

xlim([0, 2*pi])

% plot([pi, pi], [-1, 1], ':k','LineWidth',2)
% plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
% plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
% plot([0, 2*pi], [0, 0], '--k','LineWidth',0.5)
for i = 4:-2:2
    plot(outq{i}.alpha, outq{i}.current, 'LineWidth',3, 'Color', c(i, :))
end
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);