clear all

cable = dlmread('cbl_adult.dat', '\t');

%figure(1)
%hold on
%plot(cable(1:length(cable),1), cable(1:length(cable),2), 'b-')
%plot(cable(1:length(cable),1), cable(1:length(cable),3), 'g-')
%plot(cable(1:length(cable),1), cable(1:length(cable),4), 'r-')
%plot(cable(1:length(cable),1), cable(1:length(cable),5), 'c-')
%plot(cable(1:length(cable),1), cable(1:length(cable),6), 'm-')
%plot(cable(1:length(cable),1), cable(1:length(cable),7), 'y-')
%plot(cable(1:length(cable),1), cable(1:length(cable),8), 'k-')
%plot(cable(1:length(cable),1), cable(1:length(cable),9), 'b-')
%plot(cable(1:length(cable),1), cable(1:length(cable),10), 'g-')
%plot(cable(1:length(cable),1), cable(1:length(cable),11), 'r-')
%plot(cable(1:length(cable),1), cable(1:length(cable),12), 'c-')
%plot(cable(1:length(cable),1), cable(1:length(cable),13), 'm-')
%plot(cable(1:length(cable),1), cable(1:length(cable),14), 'r-')
%plot(cable(1:length(cable),1), cable(1:length(cable),15), 'c-')
%plot(cable(1:length(cable),1), cable(1:length(cable),16), 'm-')
%plot(cable(1:length(cable),1), cable(1:length(cable),17), 'y-')
%plot(cable(1:length(cable),1), cable(1:length(cable),18), 'k-')

%xlabel('t (ms)')
%ylabel('V (mV)')
%title('Voltage vs time in a Luo_Rudy cable model')

dx = 0.01; 
n = 5;

for i = 1:length(cable)
    x1(i) = n*dx;
    x2(i) = 2*n*dx;
    x3(i) = 3*n*dx;
    x4(i) = 4*n*dx;
    x5(i) = 5*n*dx;
    x6(i) = 6*n*dx;
    x7(i) = 7*n*dx;
    x8(i) = 8*n*dx;
    x9(i) = 9*n*dx;
    x10(i) = 10*n*dx;
    x11(i) = 11*n*dx;
    x12(i) = 12*n*dx;
    x13(i) = 13*n*dx;
    x14(i) = 14*n*dx;
    x15(i) = 15*n*dx;
    x16(i) = 16*n*dx;
    x17(i) = 17*n*dx;
    x18(i) = 18*n*dx;
    x19(i) = 19*n*dx;
    x20(i) = 20*n*dx;
end

figure(2)
hold on
plot3(x1, cable(1:length(cable),1), cable(1:length(cable),2), 'b-')
plot3(x2, cable(1:length(cable),1), cable(1:length(cable),3), 'g-')
plot3(x3, cable(1:length(cable),1), cable(1:length(cable),4), 'r-')
plot3(x4, cable(1:length(cable),1), cable(1:length(cable),5), 'c-')
plot3(x5, cable(1:length(cable),1), cable(1:length(cable),6), 'm-')
plot3(x6, cable(1:length(cable),1), cable(1:length(cable),7), 'y-')
plot3(x7, cable(1:length(cable),1), cable(1:length(cable),8), 'k-')
plot3(x8, cable(1:length(cable),1), cable(1:length(cable),9), 'b-')
plot3(x9, cable(1:length(cable),1), cable(1:length(cable),10), 'g-')
plot3(x10, cable(1:length(cable),1), cable(1:length(cable),11), 'r-')
plot3(x11, cable(1:length(cable),1), cable(1:length(cable),12), 'c-')
plot3(x12, cable(1:length(cable),1), cable(1:length(cable),13), 'm-')
plot3(x13, cable(1:length(cable),1), cable(1:length(cable),14), 'r-')
plot3(x14, cable(1:length(cable),1), cable(1:length(cable),15), 'c-')
plot3(x15, cable(1:length(cable),1), cable(1:length(cable),16), 'm-')
plot3(x16, cable(1:length(cable),1), cable(1:length(cable),17), 'b-')
plot3(x17, cable(1:length(cable),1), cable(1:length(cable),18), 'g-')
plot3(x18, cable(1:length(cable),1), cable(1:length(cable),19), 'r-')
plot3(x19, cable(1:length(cable),1), cable(1:length(cable),20), 'c-')
plot3(x20, cable(1:length(cable),1), cable(1:length(cable),21), 'm-')
plot3(x21, cable(1:length(cable),1), cable(1:length(cable),22), 'y-')
plot3(x22, cable(1:length(cable),1), cable(1:length(cable),23), 'k-')
plot3(x23, cable(1:length(cable),1), cable(1:length(cable),24), 'b-')
plot3(x24, cable(1:length(cable),1), cable(1:length(cable),25), 'g-')
plot3(x25, cable(1:length(cable),1), cable(1:length(cable),26), 'r-')
plot3(x26, cable(1:length(cable),1), cable(1:length(cable),27), 'c-')
plot3(x27, cable(1:length(cable),1), cable(1:length(cable),28), 'm-')
plot3(x28, cable(1:length(cable),1), cable(1:length(cable),29), 'r-')
plot3(x29, cable(1:length(cable),1), cable(1:length(cable),30), 'c-')
plot3(x30, cable(1:length(cable),1), cable(1:length(cable),31), 'm-')

xlabel('distance along cable (cm)')
ylabel('time (ms)')
zlabel('membrane voltage (mV)')
title('Voltage vs displacement vs time in a Luo-Rudy dynamic cable model')