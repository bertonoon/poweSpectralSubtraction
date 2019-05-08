clear; clc;
format compact
format long

% Parametry
koniec      = 2^18; % d³ugoœæ oryginalnego pliku audio
dl_ramka    = koniec/32; 
alfa = 0.5;
beta = 1;
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
% pause
y3 = zeros(size(y2));
for i=1:1:(koniec+dl_ramka*4)
    if i <= dl_ramka*4
        y3(i) = szum(i);
    else 
        y3(i) = y2(i-dl_ramka*4);
    end
end 
% sound(y3,Fs);


% Krok 1
ramki = 4;
S_z = zeros(dl_ramka/2,1);
N = dl_ramka;
szum_usr = zeros(dl_ramka,1);
for i=1:1:ramki
    for j=1:1:N
        szum2(j) = szum(j+(N*(i-1)));
    end
    Z = fft(szum2);
    
    for j=1:1:(N/2)
        S_z(j) = S_z(j) + (1/N)*(abs(Z(j))^2);
    end
    for j=1:1:N
        szum_usr(j) = szum2(j) +szum_usr(j);
    end
end
S_z = S_z.*1.2;
szum_usr = szum_usr./4;
% Krok 2
ramki = (length(y3)/dl_ramka);
y4 = y3(1:1:length(y3));
S_y = zeros(dl_ramka/2,1);
N = dl_ramka;
wyjsciowy = zeros(length(y4),1);

S_x = zeros(dl_ramka/2,1);
A = zeros(length(N),1);
X = zeros(length(N));
sygnal = zeros(N,1);

for i=1:1:ramki
    for j=1:1:N
        sygnal(j) = y4(j+dl_ramka*(i-1));
    end
    Y = fft(sygnal);
    
    for j=1:1:(N/2)
        S_y(j) = (1/N)*(abs(Y(j))^2);
    end 

    % Krok 3 
    for j=1:1:N/2
        if (S_y(j) - S_z(j) >= 0)
            S_x(j) = S_y(j) - S_z(i);
        else
            S_x(j) = 0;
        end
    end

    % Krok 4
    for j=1:1:N
        if j <= N/2
            A(j) = sqrt(S_x(j)/S_y(j));
        elseif j == N
            A(j) = A(1);
        else
            A(j) = A(N-j);
        end
    end
    
    for j=1:1:length(Y)
%       Metoda normalna
        X(j) = (A(j))*(Y(j));
%       Metoda rozszerzona
%         licznik = (abs(Y(j))^beta-alfa*(abs(szum_usr(j))^beta));
%         if (licznik < 0) 
%             licznik =0; 
%         end
%         X(j) = ((licznik/abs(Y(j))^beta)^(1/beta))*Y(j);
    end
    
    Xodw = ifft(X);
    for j=1:1:dl_ramka
        k = j+dl_ramka*(i-1);
        wyjsciowy(k) = Xodw(j);
    end
    i
end
wyjsciowy = real(wyjsciowy);
sound(wyjsciowy,Fs)