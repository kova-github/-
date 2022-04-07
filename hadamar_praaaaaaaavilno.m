clc;
clear;
close all;
 


SNRBitdB = -1:0.05:11;

%РС
p = zeros(1, length(SNRBitdB));% внутр РС
Pb_RS = zeros(1, length(SNRBitdB));% на бит
Pe_RS = zeros(1, length(SNRBitdB));%
%без РС
PinB = zeros(1, length(SNRBitdB)); %на бит

SNRBit = 10.^(SNRBitdB/10);


figure(1);
color = ['b','r','m','g','k','y'];
nc = 0;
%for q = [ 64,128,256,512,1024]
    for q=16

 nc = nc +1;

 N = q;
 Rin = log2(q)/q;
 for i=1:length(SNRBitdB)
    
    d = 0;
    yb = SNRBit(i);
  
    K = floor(0.99*N);
    Rout = K/N;
    
  %  for K = 1:N-1
        
        %без РС 
        ro = 0.001:0.001:1;
        tmp = ((q-1).^ro./(1+ro)).*exp(-(ro./(1+ro))*log2(q)*yb);
        Pin = min(tmp);
        
        K3 = (q/2)/(q-1);  
        PinB(i) = K3*Pin;
        
        %с РС
        tmp = ((q-1).^ro./(1+ro)).*exp(-(ro/(1+ro))*(K/N)*(log2(q)*Rout*yb));
        p = min(tmp);
                
        
        
        d = N - K + 1;
        t = (N-K)/2;
        tau = t/N;
        
        if tau < p
          Pb_RS(i) = 1;
          
        else
          Tpt = -tau*log(p) - (1-tau)*log(1-p);
          Ht = -tau*log(tau) - (1-tau)*log(1-tau);
          
          K1 = p*(1-tau)/(tau-p);
          K2 = 1/sqrt(2*pi*N*tau*(1-tau));
          
          Pe_RS = K1*K2*exp(-N*(Tpt-Ht));
          Pb_RS(i)= K3*(d/N)*Pe_RS;
        end 
   
end


%semilogy(SNRBitdB, Pb_RS,'r',SNRBitdB,PinB,'k');
%semilogy(SNRBitdB, Pb_RS,[color(nc),':'],SNRBitdB,PinB,color(nc), 'LineWidth', 2);
%semilogy(SNRBitdB,Pb_RS,[color(nc),':'], 'LineWidth', 2);
semilogy(SNRBitdB, PinB, color(nc),'LineWidth', 2);
hold on

drawnow
end
grid on;
%legend('Pb c кодом РС', 'Pb без кода РС','Location','SouthWest');
xlabel('\gamma_b, dB');
ylabel('P_b');
legend('q=64','q=128','q=256','q=512','q=1024');
%title([{'Графики вероятности ошибки на бит кода'};{' с использованием РС (...) и без использования кода РС (__)'}]);
title('Графики вероятности ошибки декодирования внутреннего кода');

axis([min(SNRBitdB), max(SNRBitdB), 1e-11, 1])


q = 128;
T = 1;% длительность или период следования сигналов 
dt = (T/q)/20;
t = 0:dt:T - dt;
 
Tc = T/q; %длительность чипа
H = hadamard(q); %Матрица Адамара
m = zeros(q, length(t)); %ортогональные огибающие
 
for i = 1:q
  tmp = repmat(H(i, :), Tc/dt, 1);
  tmp = tmp(:)';
  m(i, :) = tmp;
  k = sqrt(1/(sum(tmp.^2)*dt)); %коэффициент нормировки   
  m(i, :) = k * m(i, :);  %нормировка огибающих
end

qq = min(q,16);
for i = 1:qq
    subplot(sqrt(qq), sqrt(qq),i);
    plot(t, m(i, :), 'b');
    grid on;
    axis([-0.05, 1.05, -2, 2]);
    title(['m_{',num2str(i-1),'}(t)']);
end
 

Nerrmax = 35; % число ошибок для окончания моделирования
dBs = 1:2:13;   

Pe_pr = zeros(1, length(dBs));


SNRdBtheor = -1:0.05:11;
Pe_th = zeros(1, length(SNRdBtheor));
SNRtheor = 10.^(SNRdBtheor/10);
for l = 1:q-1
      Pe_th = Pe_th + nchoosek(q - 1, l) * (-1)^(l + 1) * 1 / (1 + l) * exp(-l/(l + 1) * SNRtheor);
end
Pe_bound = min(1, ((q - 1)/2) * exp(-SNRtheor/2));

Pb_th1 = zeros(1, length(SNRdBtheor));

Km=(q/2)/(q-1);

m = log2(q);
SNRtheorBit = SNRtheor/m;
for l = 1:q-1
      Pb_th1 = Pb_th1 + nchoosek(q - 1, l) * (-1)^(l + 1) * 1 / (1 + l) * exp(-l/(l + 1) * m*SNRtheorBit);
end
Pb_th1 = Km*Pb_th1;

Pb_bound1 = min(1,Km*((q - 1)/2) * exp(-m*SNRtheorBit/2));

for dB = 1:length(dBs)
    SNR = 10.^(dBs(dB)/10);
    Nerr = 0;
    Ntest = 0;%начальное значение счетчика числа испытаний
       
    gamma_c = SNR/q;
    sqrt2gamma_c = sqrt(2*gamma_c);
   
    %цикл моделирования при одном значении отношения сигнал/шум
    while Nerr < Nerrmax
        Ntest = Ntest + 1;
       
        %моделирование передатчика и канала
       i = floor(rand * q) + 1; %случайно выбираем i на интервале 0..q-1
       theta = 2 * pi * rand;
      
       m_i = H(i,:); 
       rc_t = cos(theta)*sqrt2gamma_c*m_i + randn(1,q);
       rs_t = sin(theta)*sqrt2gamma_c*m_i + randn(1,q);
       bpa_c=fwht(rc_t,[], 'hadamard');
       bpa_s=fwht(rs_t,[], 'hadamard');
     %  bpa_c = H*rc_t';
     %  bpa_s = H*rs_t';
       res = bpa_c .^ 2 + bpa_s .^ 2;

       
        %формирование решения
       [~, i_] = max(res);
       %фиксация результата
       if i_ ~= i
           Nerr = Nerr + 1;%увеличение счетчика числа испытаний
           clc
           disp ([dBs(dB), Ntest, Nerr, Nerrmax]);
       end
       
    end % конец цикла моделирования при одном значении сигнал/шум
    
    %вычисление экспериментальной оценки веоятности ошибки
    Pe_pr(dB) = Nerr/Ntest;
end
 
figure(2);
semilogy(SNRdBtheor, Pe_th,'g', 'LineWidth',1.5);
hold on;
%semilogy(dBs, Pe_pr, 'Color', 'm', 'LineWidth', 1.5);
semilogy(dBs, Pe_pr, 'mo', 'LineWidth', 1.5);
hold on;
semilogy(SNRdBtheor, Pe_bound, 'Color', 'k', 'LineWidth', 1.5);
grid on;
xlabel('SNR, dB');
ylabel('Pe');
legend('Pe theory', 'Pe practice');

figure(3);
SNRdBtheorBit = 10*log10(SNRtheorBit);
semilogy(SNRdBtheorBit, Pb_th1,'Color', 'b', 'LineWidth',1.5);
hold on;
%semilogy(dBs, Pe_pr, 'Color','m', 'LineWidth', 1.5);
Pb_pr = Km*Pe_pr;
SNRdBBit = dBs - 10*log10(m);
%semilogy(SNRdBBit, Pb_pr, 'mo', 'LineWidth', 1.5);
hold on;
semilogy(SNRdBtheorBit, Pb_bound1, 'Color', 'k', 'LineWidth', 1.5);
semilogy(SNRBitdB, PinB, 'Color', 'm','LineWidth', 2);
xlabel('\gamma_b, dB');
ylabel('P_b');
grid on;
title('Графики вероятности ошибки декодирования внутреннего кода, q=128');
legend('P_b точная формула', 'P_b аддитивная граница', 'P_b верхняя оценка');

