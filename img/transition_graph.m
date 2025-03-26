k = 4;
N = 3;
R = 0.6;

NT = 200;
th = linspace(0, 2*pi, NT);
c = R*cos(th);
s = R*sin(th);

figure(12)
clf
hold on
dx = 2;
dy = 2;




cc = linspecer(2);

for i=0:k-1
    for n=1:N
        up_src = [dx*i; dy*n] + R*[1; 1]*(sqrt(2)/2);
        up_dst = [dx*(i+1); dy*(n+1)] + R*[-1; -1]*(sqrt(2)/2);

        down_src = [dx*i; dy*n] +  R*[1; -1]*(sqrt(2)/2);
        down_dst = [dx*(i+1); dy*(n-1)] + R*[-1; 1]*(sqrt(2)/2);

        du = up_dst - up_src;
        dd = down_dst - down_src;

        if n < N
            quiver(up_src(1), up_src(2), du(1), du(2), 'autoscale', 1, ...
                'Color', cc(1, :), 'linewidth', 5, 'maxheadsize', 3)
        end
        if n > 1
            quiver(down_src(1), down_src(2), dd(1), dd(2), 'autoscale', 1, ...
                'Color', cc(2, :), 'linewidth', 5, 'maxheadsize', 3)
        end
    end
end


for i = 0:k
    for n = 1:N
        plot(c+dx*i, s+dy*n, 'k', 'linewidth', 3)
    end    
end

for i = 0:k
    text(dx*i-0.5, (N+0.65)*dy, sprintf('$i=%d$', i), 'fontsize', 12, 'interpreter', 'latex')
end

for n = N:-1:1
    text(-1.15*dx, dy*n, sprintf('$n=%d$', n), 'fontsize', 12, 'interpreter','latex')
end

xlim([-2*dx, (k+1)*dx])
ylim([0, (N+1)*dy])
xl = xlim;
yl = ylim;
pbaspect([diff(xl), diff(yl), 1])

axis off
