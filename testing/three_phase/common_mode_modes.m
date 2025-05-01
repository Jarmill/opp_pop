L =1;

va = kron(kron(-L:L, ones(1, 2*L+1)), ones(1, 2*L+1));
vb = kron(kron(ones(1, 2*L+1), -L:L), ones(1, 2*L+1));
vc = kron(kron(ones(1, 2*L+1), ones(1, 2*L+1)), -L:L);

V = [va; vb; vc];

Vcm3 = sum(V, 1);
ivalid = (abs(Vcm3) <= 1);

V3 = V(:, ivalid);

k2 = convhull(V3(1, :), V3(2, :), V3(3, :),'Simplify',true);

figure(1)
clf
hold on
trisurf(k2,V3(1, :),V3(2, :),V3(3, :))
scatter3(V(1, :), V(2, :), V(3, :), 100, 'k', 'filled');
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
figure(2);
clf
subplot(1, 2, 1)
plot(G)
subplot(1, 2, 2)
Gt = graph(Et);
plot(Gt)

