clear; clc;
format compact
format long

% Parametry
koniec      = 2^18; 
dl_ramka    = koniec/256; 
skal        = 10; 
ram_szum    = 30; 

% Czytanie z pliku audio
[y,Fs] = audioread('Randka w ciemno.wav');
info    = audioinfo('Randka w ciemno.wav');
y1      = y(1:1:koniec);

% Generowanie bia³ego szumu
szum = wgn(koniec,1,-40);

y2 = y1+szum;
y3 = zeros(length(y2)+dl_ramka*ram_szum,1);

for i=1:1:(koniec+dl_ramka*ram_szum)
    if i <= dl_ramka*ram_szum
        y3(i) = szum(i);
    else 
        y3(i) = y2(i-dl_ramka*ram_szum);
    end
end 

% Krok 1
ramki   = ram_szum;
S_z     = zeros(dl_ramka/2,1);
N       = dl_ramka;
szum_usr= zeros(dl_ramka,1);
szum2   = zeros(N,1);

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
ramki   = (length(y3)/dl_ramka);
y4      = y3(1:1:length(y3));
S_y     = zeros(dl_ramka/2,1);
N       = dl_ramka;
S_x     = zeros(dl_ramka/2,1);
A       = zeros(N,1);
X       = zeros(N,1);
sygnal  = zeros(N,1);
wyjsciowy = zeros(length(y4),1);


for i=1:0.5:ramki
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
        else
            A(j) = A(N-j+1);
        end
    end

    for j=1:1:length(Y)
        X(j) = (A(j))*(Y(j));
    end
    
    Xodw = real(ifft(X));
    
    for j=1:1:dl_ramka
        k = j+(N*(i-1));
        wp = 1-(abs(j-N*0.5)/(N*0.5));
        wyjsciowy(k) = wyjsciowy(k)+wp*Xodw(j);
    end
end

sound(wyjsciowy,Fs)