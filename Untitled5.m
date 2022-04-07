clc;
clear;
close all;
 

E = 1;% ������� �������� 


 
Nerrmax = 35;%��������� ����� ������
dBs = 1:14;
Pe_th = zeros(1, length(dBs));
Pe_pr = zeros(1, length(dBs));
Pe_bound = zeros(1, length(dBs));
 
figure(1);
color = ['b','r','m','k','g','y'];
nc = 0;

    for dB = 1:length(dBs)
    disp(dB);
    SNR = 10.^(dBs(dB)/10);
  
    N0 = E/SNR;

    Pe = 0;
    for l = 1:q-1
      Pe = Pe + nchoosek(q - 1, l) * (-1)^(l + 1) * 1 / (1 + l) * exp(-l/(l + 1) * E / N0);
    end
  
    Pe_th(dB) = Pe;
    Pe_bound(dB) = (q - 1) / 2 * exp(-E / (2 * N0));
  

%semilogy(dBs, Pe_th, 'LineWidth',1.5);

%semilogy(dBs, Pe_pr, 'Color', 'm', 'LineWidth', 1.5);
%hold on;
%semilogy(dBs, Pe_bound, 'Color', 'k', 'LineWidth', 1.5);
semilogy(dBs, Pe_bound(dB),color(nc), dBs, Pe_bound(dB), color(nc),'LineWidth', 2);
    
end
grid on;
%legend('Pb c ����� ��', 'Pb ��� ���� ��','Location','SouthWest');
xlabel('\gamma_b, dB');
ylabel('P_b');
legend('q=8','q=16','q=32','q=64');
%title([{'������� ����������� ������ �� ��� ����'};{' � �������������� �� (...) � ��� ������������� ���� �� (__)'}]);
title('������� ����������� ������ ������������� ����������� ���� �� ������ �������');
title([{'������� ����������� ������ �� ��� ������������� ����������� ����'};{' �� ������ �������'}]);
axis([min(SNRBitdB), max(SNRBitdB), 1e-11, 1])
 
