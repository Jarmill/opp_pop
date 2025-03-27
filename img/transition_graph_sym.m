k = 12;
N = 3;
Nc = (N+1)/2;
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

n0 = 1;

HW = -1.25*N*dy;
QW = -2.5*N*dy;

cc = linspecer(2);

%Full-Wave
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
                'Color', cc(1, :), 'linewidth', 3, 'maxheadsize', 3)
        end
        if n > 1
            quiver(down_src(1), down_src(2), dd(1), dd(2), 'autoscale', 1, ...
                'Color', cc(2, :), 'linewidth', 3, 'maxheadsize', 3)
        end
    end
end

%Half-Wave
for i=0:k/2-1
    for n=1:N
        up_src = [dx*i; dy*n + HW] + R*[1; 1]*(sqrt(2)/2);
        up_dst = [dx*(i+1); dy*(n+1) + HW] + R*[-1; -1]*(sqrt(2)/2);

        down_src = [dx*i; dy*n + HW] +  R*[1; -1]*(sqrt(2)/2);
        down_dst = [dx*(i+1); dy*(n-1) + HW] + R*[-1; 1]*(sqrt(2)/2);

        du = up_dst - up_src;
        dd = down_dst - down_src;

        if n < N
            quiver(up_src(1), up_src(2), du(1), du(2), 'autoscale', 1, ...
                'Color', cc(1, :), 'linewidth', 3, 'maxheadsize', 3)
        end
        if n > 1
            quiver(down_src(1), down_src(2), dd(1), dd(2), 'autoscale', 1, ...
                'Color', cc(2, :), 'linewidth', 3, 'maxheadsize', 3)
        end
    end
end

%Quarter-Wave
for i=0:k/4-1
    for n=1:N
        up_src = [dx*i; dy*n + QW] + R*[1; 1]*(sqrt(2)/2);
        up_dst = [dx*(i+1); dy*(n+1) + QW] + R*[-1; -1]*(sqrt(2)/2);

        down_src = [dx*i; dy*n + QW] +  R*[1; -1]*(sqrt(2)/2);
        down_dst = [dx*(i+1); dy*(n-1) + QW] + R*[-1; 1]*(sqrt(2)/2);

        du = up_dst - up_src;
        dd = down_dst - down_src;

        if n < N
            quiver(up_src(1), up_src(2), du(1), du(2), 'autoscale', 1, ...
                'Color', cc(1, :), 'linewidth', 3, 'maxheadsize', 3)
        end
        if n > 1
            quiver(down_src(1), down_src(2), dd(1), dd(2), 'autoscale', 1, ...
                'Color', cc(2, :), 'linewidth', 3, 'maxheadsize', 3)
        end
    end
end

%% vertices
%Full-wave
offset = 1.5;
text(-offset, (Nc)*dy, 'FW', 'interpreter', 'latex', 'FontSize', 16,  'HorizontalAlignment', 'right')
text(-offset, (Nc)*dy + HW, 'HW', 'interpreter', 'latex', 'FontSize', 16,  'HorizontalAlignment', 'right')
text(-offset, (Nc)*dy + QW, 'QaHW', 'interpreter', 'latex', 'FontSize', 16,  'HorizontalAlignment', 'right')


for n = 1:N
    for i = 0:k
        plot(c+dx*i, s+dy*n, 'k', 'linewidth', 3)
        if (i==0 || i==k) && n == n0
            patch(c+dx*i, s+dy*n, 'k')
        end
    end    
    for i = 0:k/2
        plot(c+dx*i, s+dy*n + HW, 'k', 'linewidth', 3)
        if (i==0 && n == n0) || (i==k/2 && n==(N-n0+1))
            patch(c+dx*i, s+dy*n+HW, 'k')
        end
    end    
    for i = 0:k/4
        plot(c+dx*i, s+dy*n+ QW, 'k', 'linewidth', 3)
        if i==0 && n==n0
            patch(c+dx*i, s+dy*n+ QW, 'k')
        end
    end    
end

% for i = 0:k
%     text(dx*i-0.5, (N+0.65)*dy, sprintf('$i=%d$', i), 'fontsize', 12, 'interpreter', 'latex')
% end
% 
% for n = N:-1:1
%     text(-1.15*dx, dy*n, sprintf('$n=%d$', n), 'fontsize', 12, 'interpreter','latex')
% end

% xlim([-2*dx, (k+1)*dx])
% ylim([0, (N+1)*dy])
xl = xlim;
yl = ylim;
pbaspect([diff(xl), diff(yl), 1])

axis off
