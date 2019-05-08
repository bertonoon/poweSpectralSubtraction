clear; clc;
format compact
format long

% Parametry
koniec      = 2^18; % d³ugoœæ oryginalnego pliku audio
dl_ramka    = koniec/512; 
skal        = 30;
skal2       = 30;
ram_szum    = 50;

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

y3 = zeros(length(y2)+dl_ramka*ram_szum,1);

for i=1:1:(koniec+dl_ramka*ram_szum)
    if i <= dl_ramka*ram_szum
        y3(i) = szum(i);
    else 
        y3(i) = y2(i-dl_ramka*ram_szum);
    end
end 

% sound(y3,Fs);
% pause

% Krok 1
ramki = ram_szum;
S_z = zeros(dl_ramka/2,1);
N = dl_ramka;
szum_usr = zeros(dl_ramka,1);
szum2 = zeros(N,1);

for i=1:1:ramki
    for j=1:1:N
        szum2(j) = szum(j+(N*(i-1)));
    end
    Z = fft(szum2);
    
    for j=1:1:(N/2)
        S_z(j) = S_z(j) + (1/N)*(abs(Z(j))^2);
    end
end
S_z = S_z./ram_szum;


% Krok 2
ramki = (length(y3)/dl_ramka);
y4 = y3(1:1:length(y3));
S_y = zeros(dl_ramka/2,1);
N = dl_ramka;
wyjsciowy = zeros(length(y4),1);
wyjsciowy2 = zeros(length(y4)-N/2,1);

S_x = zeros(dl_ramka/2,1);
A = zeros(N,1);
X = zeros(N,1);
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
        if (S_y(j) - skal*S_z(j) >= 0)
            S_x(j) = S_y(j) - skal*S_z(j);
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
%     plot(real(A));
%     pause();
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

% Overlap
for i=1:1:ramki-1
    for j=1:1:N
        sygnal(j) = y4(N/2+j+dl_ramka*(i-1));
    end
    Y = fft(sygnal);
    
    for j=1:1:(N/2)
        S_y(j) = (1/N)*(abs(Y(j))^2);
    end 

    % Krok 3 
    for j=1:1:N/2
        if (S_y(j) - skal2*S_z(j) >= 0)
            S_x(j) = S_y(j) - skal2*S_z(j);
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
%     plot(real(A));
%     pause();
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
        wyjsciowy2(k) = Xodw(j);
    end
    i
end
wyjsciowy3 = zeros(length(wyjsciowy2),1);
for i=1:1:N/2
    wyjsciowy3(i) = wyjsciowy(i);
end
k=0;

for i=1:1:length(wyjsciowy2)
    k = k+1;
    wyjsciowy3(i+N/2) = (1-abs(k-0.5*N)/(0.5*N))*wyjsciowy(i+N/2)+(1-abs(k-0.5*N)/(0.5*N))*wyjsciowy2(i);
    if k==N+1
        k=0;
    end
end

wyjsciowy3 = real(wyjsciowy3);
sound(wyjsciowy3,Fs)