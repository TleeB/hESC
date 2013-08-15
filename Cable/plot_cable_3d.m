function p=plot_cable_3d(filename,space_discretization)

% Plots action potential of a cable of cardiac cells
% [result=]plot_cable(filename, space_discretization)
% filename = name of file with cable data 
%   format of data: rows of tab-delimineted columns, where first column is the time, and
%   subsequent columns contain the voltage values at discrete spatial locations along the fiber
% space_discretization = distance (in cm) between successive spatial data points

cable = dlmread(filename, '\t');

[a,b] = size(cable);

for n = 1:b
    for i = 1:length(cable)
        x(i,n) = (n-1)*space_discretization;
    end
end

figure(2)
hold on

for n = 2:b
    plot3( x(1:length(cable),n-1), cable(1:length(cable),1), cable(1:length(cable),n))
end

view([90 70])

xlabel('X (cm)')
ylabel('Time (ms)')
zlabel('V (mV)')
