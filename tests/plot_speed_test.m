clear
clf
cpp = csvread('~/app/C++/tests/speed_test_000/speed_test.dat');
[t_mean, v,twotau] = speed_test();


figure
clf
hold on

lw = 3;
plot(twotau*1e6,real(v),'linewidth', lw);
plot(cpp(:,1)*1e6,cpp(:,2),'linewidth', lw )
yline(1)
legend('MATLAB',"C++")
box on
xlabel('2\tau (\mus)');
ylabel('coherence');
%title(['h_{18}-TEMPO in H_{2}O, no S_{z}I_{x,y}'])
fontsize = 24;
set(gca,'fontsize',fontsize);


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 12]);
print('MATLAB_vs_C++','-dpng','-r0');
set(gcf,'WindowStyle','Docked')

t_mean
% t_mean/0.00015
t_mean/0.000244

%%
[t_mean4, delta_t4] = speed_test2(4);
% cppGain4 = t_mean4/0.000138312
cppGain4 = t_mean4/0.0001387
[t_mean8, delta_t8] = speed_test2(8);
% cppGain8 = t_mean8/0.000866571
cppGain8 = t_mean8/0.0076859

[t_mean16, delta_t16] = speed_test2(8);
cppGain16 = t_mean16/0.0384927