%EN2570 Digital Signal Processing - Project

clear all;
close all;

%H.S.Watawana
%Index  - 180677E

%Deriving filter specifications from index number
Index = input('Enter your index number:');
A = mod(floor(Index/100),10);
B = mod(floor(Index/10),10);
C = mod(floor(Index),10);

Max_Ap = 0.03 + (0.01*A);            %Maximum passband ripple in dB
Min_Aa = (45 + B);                   %Minimum stopband attenuation in dB
Omega_P1 = (C*100)+300 ;             %Lower passband edge in rad/s
Omega_P2 = (C*100)+700  ;            %Upper passband edge in rad/s
Omega_S1 = (C*100)+150  ;            %Lower stopband edge in rad/s
Omega_S2 = (C*100)+800  ;            %Upper stopband edge in rad/s
Omega_Sampling = 2*((C*100)+1200);        %Sampling frequency in rad/s

%calculation of expected delta values in bands 
del_p = (10^(0.05*Max_Ap)-1) / (10^(0.05*Max_Ap)+1) ; 
del_a = (10^(-0.05*Min_Aa));

%deriving delta value for filter (ripple)
delta   = min(del_p , del_a);

%actual stop band attenuation in dB
Aa = -20*log10(delta);

%calculation of cut off frequencies and transition bandwidth
B_t = min( (Omega_S2 - Omega_P2) , (Omega_P1 - Omega_S1) );
WC2 = (Omega_P2+(B_t/2));
WC1 = (Omega_P1-(B_t/2));

%calculation of alpha
if Aa <= 21
    alpha = 0;
elseif 21<Aa && Aa <= 50
    alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21) ;
else
    alpha = 0.1102*(Aa-8.7) ;
end

%calculation of D
if Aa <= 21
    D = 0.9222;
else 
    D = (Aa - 7.95)/ 14.36;
end

%calculation of N
lower_limit=ceil((Omega_Sampling*D/B_t)+1);

if mod(lower_limit,2)==1
    N=lower_limit;
else
    N=lower_limit+1;
end

%kaiser window function
halfwidth=(N-1)/2;

%creating N wide transfer function centered around n=0
n = -halfwidth:1:halfwidth;
beta = alpha * (1 - ((2*n/(N-1)).^2)).^0.5;
I0_beta=0;
I0_alpha=0;

boundary=100;

for k=1:boundary
    I0_beta = I0_beta + ((1/factorial(k))*(beta/2).^k).^2;
    I0_alpha = I0_alpha + ((1/factorial(k))*(alpha/2)^k)^2;
end

I0_beta=I0_beta + ones(1,numel(I0_beta));
I0_alpha=I0_alpha + ones(1,numel(I0_alpha));

W=I0_beta ./ I0_alpha;

figure;
stem(n,W,'fill');
xlabel('n');
ylabel('W[n]');
title('Kaiser Windowing Function');
grid on;

%computing ideal bandpass filter impulse respose
T=2*pi/Omega_Sampling;
width=(N-1)/2;
n_neg=-width:1:-1;
hI_neg= ((1./n_neg)*(1/pi)).*(sin(WC2*n_neg*T)-sin(WC1*n_neg*T));
hI_0=2*(WC2-WC1)/Omega_Sampling;
n_pos=1:1:width;
hI_pos=((1./n_pos)*(1/pi)).*(sin(WC2*n_pos*T)-sin(WC1*n_pos*T));


hI=[hI_neg,hI_0,hI_pos];
n=[n_neg,0,n_pos];

figure;
stem(n,hI,'fill');
xlabel('n');
ylabel('hI[n]');
title('Ideal Impulse Response of a Bandpass Filter');
grid on;


%Computing impulse response for FIR bandpass filter
h=hI.*W; 

figure;
stem(n,hI,'fill');
xlabel('n');
ylabel('h[n]');
title('Impulse Response of a FIR Non-Causal Bandpass Filter');
grid on;

shift_n=[0:1:N-1];
figure;
stem(shift_n,h,'fill');
xlabel('n');
ylabel('hI[n]');
title('Impulse Response of a FIR Causal Bandpass Filte1r');
grid on;

%Magnitude response of the Bandpass Filter
fvtool(h);
freqz(h);

%Output of filter to a sudden excitation
Omega1=Omega_S1/2;
Omega2=(Omega_P1+Omega_P2)/2;
Omega3=(Omega_S2+Omega_Sampling)/2;

n=0:1:300;
m=length(n);
x=sin(Omega1*n*T)+sin(Omega2*n*T)+sin(Omega3*n*T);
L=length(h);

X=[x,zeros(1,L)];
H=[h,zeros(1,m)];

for i=1:L+m-1
    Y(i)=0;
    for j=1:m
        if (i-j+1>0)
            Y(i)=Y(i)+X(j)*H(i-j+1);
        end
    end
    
end

figure;
stem(X,'fill');
xlabel('n');
ylabel('X[n]');
title('Excitation');
grid on;

figure;
stem(Y,'fill');
xlabel('n');
ylabel('Y[n]');
title('Output Signal');
grid on;

H_I=[hI,zeros(1,m)];

for i=1:L+m-1
    Y(i)=0;
    for j=1:m
        if (i-j+1>0)
            Y(i)=Y(i)+X(j)*H_I(i-j+1);
        end
    end
    
end

figure;
stem(Y,'fill');
xlabel('n');
ylabel('Y[n]');
title('Output Signal from Ideal Filter');
grid on;


%plotting DFT of the excitation signal
m = numel(n);
nf = 2^nextpow2(m) ; % Next power of 2 from length of y
Y = fft(x ,nf)/m;
f = (Omega_Sampling) /2* linspace (0 ,1 ,nf/2+1) ;
figure ;
subplot (3 ,1 ,1) ;
plot( f ,2*abs (Y( 1 :nf/2+1) ) ) ;
title('DFT of the excitation signal')
xlabel ('frequency in radians per second');
ylabel ( '|x(f)|' ) ; 
grid on;

%Plot the DFT of filtered signal

x_f = conv ( x , h ,'same' ) ;
m = numel ( x_f ) ;
nf = 2^nextpow2(m) ; % Next power of 2 from length of y
Y = fft( x_f ,nf) /m ;
f = (Omega_Sampling) /2* linspace (0 ,1 ,nf/2+1) ;
subplot (3 ,1 ,2) ;
plot( f ,2* abs (Y( 1 :nf/2+1) ) ) ;
title('DFT of the filtered signal' );
xlabel ('Frequency ( rad/ s )' );
ylabel (' |X( f ) | ') ; 
grid on;



%Plot the excitation in time domain

figure ;
stem(n, x ,'r','fill' ) ;
title('Excitation in time domain')
xlabel ('n')
ylabel ('x[n]') ; 
grid on;

%Plot the Filtered signal in the time domain

nl = [0 : 1 :numel( x_f )-1];

figure ;
stem( nl , x_f , 'r', 'fill' ) ;
title( 'filtered sgnal in time domain');
xlabel ('n' )
ylabel ('x[n]') ; 
grid on;

%Plot the Ideally Filtered s ignal in the time domain
x_i =conv ( x , hI ,'same' );
figure ;
stem(n, x_i ,'r','fill' ) ;
title('signal filtered using the ideal filter in time domain');
xlabel ('n');
ylabel ('x[n]') ; 
grid on;







    
    
















