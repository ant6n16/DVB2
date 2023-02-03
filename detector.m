function bits_rx = detector (datos_rx,CONSTEL)

% Detector generico para QPSK, 16QAM y 64QAM
I = real(datos_rx); Q = imag(datos_rx);

switch CONSTEL
    
    case 'BPSK'
        b0 = I; b0 = b0<0;
        bits_rx = b0;

    
    case 'QPSK'
        
        b0 = I; b1 = Q; b0 = b0<0; b1 = b1<0;
        bits_rx = [b1, b0];
        
  
    case '16QAM'
        
        escala = 2/sqrt(10);
        b0 = I; b1 = Q; b2 = abs(I) - escala; b3 = abs(Q) - escala;
        bits_rx = [b0,b1,b2,b3];
        bits_rx = bits_rx<0;
        
    case '64QAM'
        
        escala1 = 4/sqrt(42); escala2 = 2/sqrt(42);
        b0 = I; b1 = Q; b2 = abs(I) - escala1; b3 = abs(Q) - escala1;
        b4 = abs(b2) - escala2; b5 = abs(b3) - escala2;
        bits_rx = [b0,b1,b2,b3,b4,b5];
        bits_rx = bits_rx<0;
        
end


end

