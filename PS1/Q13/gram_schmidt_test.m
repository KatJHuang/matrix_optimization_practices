clc,clear,close all
v = [3 2;1 2];
dx = 1;

h=figure;hold on
for ii=1:2
    plotVector(v(:,ii),h,'r--o');
end
W = gram_schmidt(v, dx);
for ii=1:2
    plotVector(W(:,ii),h,'b');
end
axis equal
grid on;box on;

    function plotVector(V,h,C)
        % a function to plot a vector V on a figure given by handle h
        % C specifies the line properties
        figure(h);
        V = V(:);
        V = V / norm(V);
        plot([0,V(1,1)],[0,V(2,1)],C,'linewidth',1.5);
    end