load('pagerank_adj.mat');
A = get_freq(J);

[u1, error] = power_iterate(A);

plot(error)