d = 7;
N = 5;
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


start_level =3;

cc = linspecer(2);

for i=0:d-1
    for n=1:N
        % if start_level==0 || i >= (N-1)/2 || 
        if true
            up_src = [dx*i; dy*n] + R*[1; 1]*(sqrt(2)/2);
            up_dst = [dx*(i+1); dy*(n+1)] + R*[-1; -1]*(sqrt(2)/2);
    
            down_src = [dx*i; dy*n] +  R*[1; -1]*(sqrt(2)/2);
            down_dst = [dx*(i+1); dy*(n-1)] + R*[-1; 1]*(sqrt(2)/2);
    
            du = up_dst - up_src;
            dd = down_dst - down_src;
    
            if n < N && (mod(i+n, 2)==1) && (start_level==0 || i >= (N-1)/2 || abs(n-start_level)<=i)
                quiver(up_src(1), up_src(2), du(1), du(2), 'autoscale', 1, ...
                    'Color', cc(1, :), 'linewidth', 5, 'maxheadsize', 3)
            end
            if n > 1 &&  (mod(i+n, 2)==1) &&(start_level==0 || i >= (N-1)/2 || abs(n-start_level)<=i)
                quiver(down_src(1), down_src(2), dd(1), dd(2), 'autoscale', 1, ...
                    'Color', cc(2, :), 'linewidth', 5, 'maxheadsize', 3)
            end
        end
    end
end


for i = 0:d
    for n = 1:N
        if  (mod(i+n, 2)==1) && (start_level==0 || i >= (N-1)/2 || abs(n-start_level)<=i)
            plot(c+dx*i, s+dy*n, 'k', 'linewidth', 3)
        end
    end    
end

for i = 0:d
    text(dx*i-0.5, (N+0.65)*dy, sprintf('$i=%d$', i), 'fontsize', 12, 'interpreter', 'latex')
end

for n = N:-1:1
    text(-1.15*dx, dy*n, sprintf('$n=%d$', n), 'fontsize', 12, 'interpreter','latex')
end

xlim([-2*dx, (d+1)*dx])
ylim([0, (N+1)*dy])
xl = xlim;
yl = ylim;
pbaspect([diff(xl), diff(yl), 1])

axis off
