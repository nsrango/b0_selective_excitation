%%
clear
load('siemens_standard_rf90');
rf = siemens_standard_rf90' /sum(siemens_standard_rf90) * (pi/2);
x = -2 *3.5:0.01:2 *3.5;
ab  = abr(rf, x);
profile = abs(ab2ex(ab));
plot(toHz(x, 3),profile)
title('HASTE 90 degree pulse on Siemens scanner')
xlabel('Hz')

%% SlR BW illustration 
clear
%rf0 = dzrf(500, 4.5, 'st', 'ms');
%rf0 = (pi) * rf0 / sum(rf0);
%rf1 = (pi) * msinc(500, 4.5/4);
rf2 = dzrf(500, 4.5, 'ex', 'pm',0.0006, 0.012);
rf3 = dzrf(500, 4.5, 'ex', 'pm',0.01, 0.012);
%rf2 = (pi/18) * rf2 /sum(rf2);
%rf2 = dzrf(500, 4.5, 'inv', 'pm',0.0006, 0.012);
%%
%rf2 = (pi/2) * rf2 /sum(rf2);
x = -4 *3.5:0.01:4 *3.5;
ab2  = abr(rf2, x);
figure(3)
plot(toHz(x, 3),abs(ab2ex(ab2)))
hold on 

ab3  = abr(rf3, x);
%figure(4)
aaa = abs(ab2ex(ab3));
plot(toHz(x, 3),abs(ab2ex(ab3)))
title('Illustration of SLR')
legend('B = 1.5kHz' ,'B = 5kHz, shifted')
hold off;
%% Fourier BW illustration
clear
rf2 = (pi/2) * msinc(500, 4.5/4);
rf3 = (pi/2) * msinc(500, 15/4);
x = -4 *3.5:0.01:4 *3.5;
ab2  = abr(rf2, x);
figure(3)
plot(toHz(x, 3),abs(ab2ex(ab2)))
hold on 

x1 = -6 *3.5:0.01:6 *3.5;
ab3  = abr(rf3, x1);
%figure(4)
aaa = abs(ab2ex(ab3));
plot(toHz(x, 3),aaa(201:3001))
title('Illustration of Fourier')
legend('B = 1.5kHz' ,'B = 5kHz, shifted')


%%
x = -2 *3.5:0.01:2 *3.5;
ab0  = abr(rf0, x);
ab1  = abr(rf1, x);
ab2  = abr(rf2, x);
figure(1)
plot(toHz(x, 3),real(ab2inv(ab0)),toHz(x, 3),imag(ab2inv(ab0)))
title('slr small tip, optimized by msinc, 180 degree profile')
saveas(gcf,'results\41111.png')
figure(2)
plot(toHz(x, 3),real(ab2inv(ab1)),toHz(x, 3),imag(ab2inv(ab1)))
title('Fourier small tip, 180 degree profile')
saveas(gcf,'results\51111.png')
figure(3)
plot(toHz(x, 3),real(ab2inv(ab2)),toHz(x, 3),imag(ab2inv(ab2)))
title('slr excitation, optimized by PM, 180 degree profile')
saveas(gcf,'results\61111.png')

%%
clear;
TB = 4.5;
duration = 3;
%rf2 = dzrf(500, TB, 'ex', 'pm',0.0008, 0.014);
rf2 = dzrf(500, TB, 'ex', 'pm',0.02, 0.03);

x = -2 * 3.5: 0.01 :2 *3.5;
%x = -10:0.01:10; %large TB ,e.g. 15
ab2  = abr(rf2, x);
figure(1)
profile = abs(ab2ex(ab2));
save(['/home/molin/shim/results/newripplesTB_' num2str(TB), '_t_' num2str(duration) '.mat'], 'profile')
plot(toHz(x, duration),profile)
title(['slr excitation, optimized by PM, TB = ' num2str(TB), ', t = ' num2str(duration) 'ms'])
saveas(gcf,['results\' num2str(TB) '_' num2str(duration) '.png'])

%% sub-nyquist
TB = 10;
dur = 0.5;
sampled = 30;
rf_sub = dzrf(60, TB, 'ex', 'pm',0.01, 0.01);
rf0 = dzrf(sampled, TB, 'st', 'ms',0.01,0.01);
rf0 = (pi/2) * rf0 / sum(rf0);

x = -1000000 : 10 : 1000000;
ab1  = abr(rf_sub, x);
ab2 = abr(rf0,x);
figure(1)
profile = abs(ab2ex(ab1));
plot(toHz(x, dur),profile)
title(['SLR of subsampling ' num2str(sampled)])
xlabel('Hz')
figure(2)

profile1 = abs(ab2ex(ab2));
plot(toHz(x, dur),profile1)
ylim([0,1])
title(['Small tip angle of subsampling ' num2str(sampled)])
xlabel('Hz')
%%
function H = toHz(x, duration)
    H = 1000 * x/duration; 
end
