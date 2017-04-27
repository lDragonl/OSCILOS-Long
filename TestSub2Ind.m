N = 100000;
for ii = 1:N;
tic
    sub2ind(size(dataGhost), 17, 42);
a(ii) = toc;
end

figure(1)
clf
hist(a,round(sqrt(N)))


M = size(dataGhost, 1);
for ii = 1:N;
tic
    idx = M * (42 - 1) + 17;
b(ii) = toc;
end

hold on
hist(b,round(sqrt(N)))