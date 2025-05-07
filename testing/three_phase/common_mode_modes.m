L =1;

va = kron(kron(-L:L, ones(1, 2*L+1)), ones(1, 2*L+1));
vb = kron(kron(ones(1, 2*L+1), -L:L), ones(1, 2*L+1));
vc = kron(kron(ones(1, 2*L+1), ones(1, 2*L+1)), -L:L);

V = [va; vb; vc];

Vcm3 = sum(V, 1)/3;
ivalid = (abs(Vcm3) <= 1/3);

V3 = V(:, ivalid);

sl = {'-', '0', '+'};
nl = cell(size(V3, 2), 1);
for i = 1:length(V3)
    si = strcat(sl(V3(1, i)+2), sl(V3(2, i)+2), sl(V3(3, i)+2));
    nl{i} = si{1};
end

k2 = convhull(V3(1, :), V3(2, :), V3(3, :),'Simplify',true);

figure(1)
clf
tiledlayout(1, 2)
nexttile
hold on
trisurf(k2,V3(1, :),V3(2, :),V3(3, :), 'EdgeColor', 'k', 'FaceColor', 0.6*[1, 1, 1],...
    'FaceAlpha', 0.8)
scatter3(V(1, :), V(2, :), V(3, :), 100, 'k', 'filled');
xlabel('$u_a$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$u_b$', 'interpreter', 'latex', 'FontSize', 14)
zlabel('$u_c$', 'interpreter', 'latex', 'FontSize', 14)
axis equal
view(3)

N3 = size(V3, 2);
Eh = zeros(N3);
Ep = zeros(N3);
for i = 1:N3
    for j = 1:(i-1)
        dv = abs(V3(:, i) - V3(:, j));
        if max(dv) <= 1
            Eh(i, j) = 1;
        end
        if nnz(dv)==1
            Ep(i, j) = 1;
        end
    end
end
E = Eh + Eh';
Et = Ep + Ep';

E1 = E(V3(1, :)==-1, V3(1, :)==-1);
E2 = E(V3(1, :)==0, V3(1, :)==0);
E3 = E(V3(1, :)==1, V3(1, :)==1);


G = graph(E);
% figure(2);
% clf
nexttile
% subplot(1, 2, 1)
% plot(G)
% subplot(1, 2, 2)
Gt = graph(Et);
h = plot(Gt, 'nodelabel', nl);
h.NodeFontSize = 12;
axis square
axis off

