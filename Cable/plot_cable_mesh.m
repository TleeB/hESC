function p=plot_cable_mesh(filename,step)

% Plots action potential of a cable of cardiac cells
% [result=]plot_cable(filename,step)
% filename = name of file with cable data
% step = number of data data points to skip

if nargin==1
step=1;
end

data=dlmread(filename, '\t');
s=size(data);

mesh(data(1:step:s(1,1),2:s(1,2)));
view(90,75);

grid off;
colormap(gray);
p=0;
