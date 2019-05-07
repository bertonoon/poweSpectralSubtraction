clear; clc;
format compact
format long

% Parametry
koniec      = 2^18; % d³ugoœæ oryginalnego pliku audio
dl_ramka    = koniec/32; 

% Czytanie z pliku audio
[y,Fs] = audioread('Randka w ciemno.wav');
info = audioinfo('Randka w ciemno.wav');
t = 0:seconds(1/Fs):seconds(info.Duration);
t = t(1:end);

t1 = t(1:1:koniec);
y1 = y(1:1:koniec);

% sound(y1,Fs)

% Generowanie bia³ego szumu
szum = wgn(koniec,1,-30);
% sound(szum,Fs);
% plot(1:1:Fs,szum);

y2 = y1+szum;
% sound(y2,Fs)
y3 = zeros(size(y2));

for i=1:1:(koniec+dl_ramka*4)
    if i <= dl_ramka*4
        y3(i) = szum(i);
    else 
        y3(i) = y2(i-dl_ramka*4);
    end
end
% sound(y3,Fs);

y4 = y3(Fs:1:koniec+dl_ramka*4);
Y = fft(y4);

% Krok 1
Z = fft(szum);
ramki = 4;
S_z = zeros(dl_ramka/2,1);
N = dl_ramka;
S_z_suma = 0;
for i=1:1:ramki
    k = 1;
    for j=1:1:N
        szum2(j) = szum(j+i*(i-1));
    end
    Z = fft(szum2);
    
    for j=1:1:(N/2)
        S_z(j) = S_z(j) + abs(Z(j))^2;
    end 
end
S_z = S_z./ramki;

% Krok 2
Y = fft(y4);
ramki = koniec/dl_ramka;
% S_z = zeros(dl_ramka/2,1);
% N = dl_ramka;
% 
% for i=1:1:ramki
%     k = 1;
%     for j=(i*N-(N-1)):1:(i*N-N/2)
%         S_z(k) = S_z(k) + abs(Z(j))^2;
%         k = k+1;
%     end
% end
S_z = S_z./N;


