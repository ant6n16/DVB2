clear all;
close all; clc
%% Configuracion del sistema OFDM

MODO = '2K'; % Elegir 2K u 8K
CP = 1/4; % Elegir 1/4, 1/8, 1/16, 1/32
CHANNEL = 'P1'; % Elegir recepcion portable: P1 o recepcion fija: F1

ESTIMADOR = 'MATLAB'; % Elegir MATLAB si se quiere estimar el canala en 
% MATLAB o VHDL si se quiere importar la estimacion realizada en VHDL

switch MODO
    case '2K'
        NFFT = 2048; % Numero de portadoras de la OFDM
        NCARRIER = 1705; % Numero de portadoras utilizadas
        Tu = 224e-6; % Tiempo entre portadoras
        
       
    case '8K'
        NFFT = 8192; % Numero de portadoras de la OFDM
        NCARRIER = 6817;% Numero de portadoras utilizadas 
        Tu = 896e-6; %Tiempo entre portadoras
end

NPILOTOS = ceil(NCARRIER/12);
NDATA = NCARRIER - NPILOTOS;
NCP = NFFT*CP; % Numero de muestras del prefijo ciclico
                     
NUM_SYMB = 10;       % Numero de simbols a transmitir - 10 simbolos OFDM 
                     % con su prefijo ciclico
SEED=100;            % Semilla para el generador de numeros aleatorios

CONSTEL = 'QPSK';    % Constelacion utilizada BPSK o QPSK
SNR=100;             % SNR en dB

tic

% Inicializamos el generador de numeros aleatorios con la semilla -> Puede
% ser uniforme, gaussiana... al usar el generador siempre saldran los
% mismos valores aleatorios al usar la misma semilla de ruido
rng(SEED);

% Definicion de la constelacion
switch CONSTEL
    case 'BPSK'
        M = 1;    % Numero de bits
        C = [1 -1];
    case 'QPSK'
        % Entre raiz de dos para normalizar energia
        C = [1+1i 1-1i -1+1i -1-1i]/sqrt(2);
        M = 2;    % Numero de bits - 4 puntos transmitidos en 2 bits
        
    case '16QAM'
        M = 4;          
        C =  [3+3i -3+3i 3-3i -3-3i 1+3i -1+3i 1-3i -1-3i 3+1i -3+1i...
              3-1i -3-1i 1+1i -1+1i 1-1i -1-1i]/sqrt(10);

        % mean(abs(C).^2) -> Para comprobar el factor de normalizacion
        
    case '64QAM'
        M = 6;          
        C =  [7+7i -7+7i 7-7i -7-7i 1+7i -1+7i 1-7i -1-7i 7+1i -7+1i...
              7-1i -7-1i 1+1i -1+1i 1-1i -1-1i 5+7i -5+7i 5-7i -5-7i...
              3+7i -3+7i 3-7i -3-7i 5+1i -5+1i 5-1i -5-1i 3+1i -3+1i...
              3-1i -3-1i 7+5i -7+5i 7-5i -7-5i 1+5i -1+5i 1-5i -1-5i...
              7+3i -7+3i 7-3i -7-3i 1+3i -1+3i 1-3i -1-3i 5+5i -5+5i...
              5-5i -5-5i 3+5i -3+5i 3-5i -3-5i 5+3i -5+3i 5-3i -5-3i...
              3+3i -3+3i 3-3i -3-3i]/sqrt(42);
end

scatterplot(C);
grid
title('Constelacion')

%% Transmisor
% Generacion de los bits a transmitir
numbits = NUM_SYMB*NDATA*M; % Numero simbolos por Numero portadoras 
                            % ocupadas por numero bits
                            
bits_tx = rand(numbits, 1)>0.5; % numbits x 1 -> Lo pasa a bits (0 o 1)

% Bits to symbols
aux  = reshape(bits_tx, M, []).'; % numbits/M x M
symb = zeros(size(aux, 1),1);     % numbits/M x 1
for k=1:M
    symb = symb + (2^(k-1))*aux(:,k); % primera columna = lsb 
end

% Mapper
const_points = C(symb+1); % numbits/M x 1 -> le sumo uno porque en MATLAB 
                          % los indices empiezan en uno

scatterplot(const_points);
grid
title('Constelacion transmitida')

%% PRBS -> Pilotos:
vector = ones(1,11);  % Vector de unos inicial
out = zeros(1,NCARRIER); % Secuencia de salida, un bit por portadora
pilots_pos = 1:12:NCARRIER; % Donde van los pilotos dentro del bloque

% Un bit por piloto pero por estandar se genera una secuencia con un bit
% para cada portadora. De esta secuencia cogere uno de cada doce para
% introducir los pilotos

for ii=1:NFFT   
    out(ii) = vector(11);
    out_XOR = bitxor(vector(9),vector(11));
    vector = [out_XOR vector(1:10)];
end

% Modular los pilotos: BPSK
out_interes = out(1:12:NCARRIER);
Pilotos_mod = 4/3*2*(1/2-out_interes);
Pilotos_bloque = repmat(Pilotos_mod,NUM_SYMB,1); % Para estimar el canal 
% en cada simbolo  

% TPS:
tps = zeros(NCARRIER,1);
tps(1:12:NCARRIER) = Pilotos_mod;

% Repetimos para cada simbolo:
bloque = repmat(tps,1,NUM_SYMB);

% Relleno el bloque con los bits a transmitir
datos = reshape(const_points, NDATA, NUM_SYMB);
bloque(bloque==0) = datos(:);

%% Simbolos OFDM en frecuencia (rejilla tiempo frecuencia)
ofdm_freq = zeros(NFFT, NUM_SYMB); % NFFT x NUM_SYMB = 128 x 10
                                   % filas frecuencias, columnas simbolos
                                   % OFDM
                                                                   
% Esto es para construir como el espectro. Las primeras 16 filas son ceros 
% ya que esas subportadoras no se usan. Empiezo en la 17. Los puntos de la
% constelacion estan en un vector, hay que ponerlos de forma adecuada. Se 
% ve bien en la representacion que hay justo aqui abajo:                                    
% ofdm_freq(ceil((NFFT-NCARRIER)/2)+(1:NCARRIER),:) = reshape(const_points, NCARRIER, NUM_SYMB);
ofdm_freq(ceil((NFFT-NCARRIER)/2)+(1:NCARRIER),:) = reshape(bloque, NCARRIER, NUM_SYMB);

figure
stem(abs(ofdm_freq(:,1))); % Pintamos un unico simbolo
grid
xlabel('Portadoras OFDM');
ylabel('Amplitud');
title('Espectro OFDM')

% ifftshift permite pasar de una representacion del espectro con el f=0 en 
% el centro a una representacion con f=0 a la izquierda.
% Importante el 1 para hacer la transformacion en la dimension correcta
% (por columnas)
ofdm_freq=ifftshift(ofdm_freq, 1); % NFFT x NUM_SYMB

% Modulacion OFDM -> Paso de frecuencia a tiempo
% Importante el 1 para hacer la transformacion en la dimension correcta
ofdm_time = ifft(ofdm_freq, NFFT, 1); % NFFT x NUM_SYMB

% Prefijo ciclico de forma matricial
ofdm_time = [ofdm_time(end-(NCP-1):end, :); ofdm_time];

% Salida secuencial (el : lee por columnas)
tx = ofdm_time(:); % (NFFT+NCP)??NUM_SYMB x 1

figure
plot(real(tx), 'b-');
hold on
plot(imag(tx), 'r-');
xlabel('Muestras temporales');
ylabel('Amplitud');
legend('real', 'imag');
grid
title('Senal OFDM en el tiempo')

% Espectro de la senal transmitida -- Estimador densidad espectral potencia
figure
pwelch(tx);

%% CANAL P1

% VALORES DE LA TABLA: Cada entrada (fila) es una delta (20 deltas)
rho = [0.057662 0.176809 0.407163 0.303585 0.258782 ...
    0.061831 0.150340 0.051534 0.185074 0.400967 ...
    0.295723 0.350825 0.262909 0.225894 0.170996 ...
    0.149723 0.240140 0.116587 0.221155 0.259730];

tau = (1e-6).*[1.003019 5.422091 0.518650 2.751772 0.602895 ...
               1.016585 0.143556 0.153832 3.324866 1.935570 ...
               0.429948 3.228872 0.848831 0.073883 0.203952 ...
               0.194207 0.924450 1.381320 0.640512 1.368671];
 
theta = [4.855121 3.419109 5.864470 2.215894 3.758058 ...
        5.430202 3.952093 1.093586 5.775198 0.154459 ...
        5.928383 3.053023 0.628578 2.128544 1.099463 ...
        3.462951 3.664773 2.833799 3.334290 0.393889 ];
    
% Valor de k:
k = 1/sqrt(sum(rho.^2));

% Para 2K y ancho de banda de 8MHZ la separacion entre portadoras es de:
sep_portadoras = 1/(Tu); % Aprox 4khz, diferente en 8k
Tm = (Tu)/NFFT; % 2048 portadoras en el modo 2K
fm = NFFT/(Tu); 

% Vector de frecuencias:
f = ((1:NFFT)-(NFFT/2+1))*sep_portadoras; % Lo centro con el -1025.
w = 2*pi*f; % Esto es la w discretizada

switch CHANNEL
    
    case 'P1'

        HP1 = 0; % Canal al principio a cero

        % Para cada una de las 20 deltas usamos la expresion del estandar,
        % cambiando la delta por una exponencial
        for ff = 1:length(rho)    
            y = k*rho(ff)*exp(-1j*theta(ff))*exp(-1j*w*tau(ff));
            HP1 = HP1 + y;  
        end

        % Representacion del canal 
        figure, plot(f,20*log10(abs(HP1)))
        xlabel('f(MHz)'); ylabel('|H| (dB)'); grid on
        axis([min(f) max(f) -40 10]), title('Canal P1')
        
        
    case 'F1'
        
        rho_0 = sqrt(10)*1/k; 
        den = sqrt(rho_0.^2 + sum(rho.^2));
        HF1 = rho_0; % Canal al principio, rayo Rice

        % Para cada una de las 20 deltas usamos la expresion del estandar,
        % cambiando la delta por una exponencial
        for ff = 1:length(rho)   
            y = (exp(-1j*theta(ff))*exp(-1j*w*tau(ff)))/den;
            HF1 = HF1 + y;  
        end

        % Representacion del canal 
        figure, plot(f,20*log10(abs(HF1)))
        xlabel('f(MHz)'); ylabel('|H| (dB)'); grid on, title('Canal F1')
        HP1 = HF1;
             
end

%% Pasar canal al dominio del tiempo y hacer convolucion con la senal tx:
HP1_temp = ifft(fftshift(HP1),NFFT);
tx_P1 = conv(tx,HP1_temp);
tx_P1 = tx_P1(1:end-(length(tx_P1)-length(tx)));

%% Canal AWGN

% Ruido
Ps = mean(tx.*conj(tx)); % Potencia de senal
nsr = 10^(-SNR/10);      % Pn/Ps -> Potencia de ruido entre la de senal.

noise = (randn(size(tx))+1i*randn(size(tx))) / sqrt(2); % Ruido complejo de potencia 1
noise = sqrt(Ps*nsr).*noise; % Ruido complejo de potencia Ps*snr
% Alternativa a las dos lineas anteriores:
% noise = wgn(size(tx,1), 1, Ps*nsr, 'complex');

rx = tx_P1 + noise;

%% Receptor
% TO-DO: Receptor (operaciones inversas al transmisor)

% Paso el vector a matriz
rx_ofdm = reshape(rx,NFFT+NCP,NUM_SYMB);

% Quitar prefijo ciclico:
ofdm_rec = (rx_ofdm(NCP+1:end,:));

% Paso a frecuencia:
ofdm_rec_f = fft(ofdm_rec, NFFT, 1);

% Modificacion visualizacion del espectro con fftshift:
ofdm_rec_f = fftshift(ofdm_rec_f, 1);

% Representamos 1 simbolo (t=1) en el dominio de la frecuencia
figure
stem(abs(ofdm_rec_f(:,1)));
grid
xlabel('Portadoras OFDM');
ylabel('Amplitud');
title('Espectro OFDM Recibido')

% Quitar los ceros anadidos arriba y abajo:
constelacion_rx = ofdm_rec_f(ceil((NFFT-NCARRIER)/2)+(1:NCARRIER),:);

%% Interpolador - Estimamos la variacion lineal del canal entre los pilotos

switch ESTIMADOR
    
    case 'MATLAB'   
        H_estimada = zeros(size(bloque)); % Canal estimado

        Pilotos_rx = constelacion_rx(1:12:NCARRIER,:); % Extraer los pilotos
        H_estimada(pilots_pos,:) = Pilotos_rx./Pilotos_bloque'; % Pilotos ecualizados

        for pp = 1:NPILOTOS-1
            for x = 1:12
                i = x + pilots_pos(pp);
                H_estimada(i,:) = ...
                    H_estimada(pilots_pos(pp),:)*((12-x)/12)+ ...
                    (H_estimada(pilots_pos(pp+1),:)*(x/12));  
            end
        end
        
    case 'VHDL'
        Estim_VHDL = load('Verificacion/FICHEROS_MAT/Canal_Estimado_VHDL');
        H_estimada = [Estim_VHDL.H_est; 0+0j]; % Para tener 1705
end

%% Representacion del canal estimado sobre el real

% Representacion del canal real
HP1_est = zeros(size(HP1));
Dif = ceil((NFFT-NCARRIER)/2);
HP1_est(1,Dif:(NFFT-Dif)) = H_estimada(:,1);
figure, plot(f,20*log10(abs(HP1))), hold on, grid on
plot(f,20*log10(abs(HP1_est)),'r'), axis([-5e6 5e6 -40 10]),
title('Canal real y canal estimado '),legend('Canal real','Canal estimado')

%% Ecualizacion
datos_rx = constelacion_rx./H_estimada;
datos_rx(pilots_pos,:) = []; % Quito los pilotos para tener los datos solo
datos_rx = datos_rx(:); % Los pongo como vector

%% Representacion constelacion recibida

scatterplot(datos_rx);
grid
title('Constelacion recibida')
axis([-2 2 -2 2]), hold on, plot(C,'C*','MarkerSize',10)
legend('Constelacion recibida', 'Constelacion enviada')

%% Demap

bits_rx = zeros(length(datos_rx)*M,1);
bits_rx = reshape(bits_rx,size(aux));
for ii=1:length(datos_rx)          
    bits_rx(ii,:) = detector(datos_rx(ii),CONSTEL);                  
end
                  
bits_rx = reshape(bits_rx.',length(datos_rx)*M,1);

BER = mean(xor(bits_rx, bits_tx));
fprintf(1, 'BER = %f\n', BER);

% El programa para aqui, si se quiere generar los ficheros y salvar las
% variables correspondientes para verificacion ejecutar la celda posterior
return
