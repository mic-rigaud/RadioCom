clear all, close all, clc,

% Nombre d'états de la QAM.
M = 4;
% Nombre de porteuses dans le symbole OFDM
Nb = 4;
% Nombre de symboles OFDM dans la simulation
NbSym = 100; %100 symboles 
% Tirage aléatoire d'entiers allant de 0 à M-1
R = randint(Nb*NbSym,1,M);
% Mise en constellation QAM.
x1 = pskmod(R, M);

% insertion pilote 1+j tous les 3 symboles OFDM
pilote=[1+1j;1+1j;1+1j;1+1j];
piloteIFFT=ifft(pilote);
X=[];
NbSymNew=NbSym+NbSym/3;
bou=1;
for i = 1:NbSym+NbSym/3
if mod(i,3)==0
    X=[X;pilote];
else
   X=[X;x1((bou-1)*4+1:bou*4)];
   bou=bou+1;
end
    
end
scatterplot(X);title('Constelations du signal à l emission après modulation pi/4-QPSK');figure;

% Création signal OFDM
x = zeros(size(X));
for ind = 1:NbSymNew
% calcul ième symbole OFDM
symbole=ifft(X((ind-1)*Nb+1:ind*Nb));
% sauvegarde du symbole ind dans x
x((ind-1)*Nb+1:ind*Nb) = symbole; % bien comprendre que l'info x(1) est sur la fréquence de la sous-porteuse 1, x(2) sur la 2,..., x(5) sur la 1,...
end
subplot(2,1,1); plot(real(x))
title('partie réelle de x')
subplot(2,1,2); plot(imag(x))
title('partie imaginaire de x')

%% Canal de propagation
figure;
% passage dans le canal multi-trajet
% On a mit des fréquences comme pour le canal 1 du wifi
% Norme 802.11a : espacement entre sous-porteuses 0.3125MHz
fcentre=2.412e9;
espacementSP=0.3125e6;
f=[fcentre,fcentre+espacementSP,fcentre+2*espacementSP,fcentre+3*espacementSP]; 
h=ones(1,4)+(0.4+j*0.2).*(1./f); % réponse fréquenciel d'un canal écho
subplot(2,1,1)
plot(f,abs(h)) % réponse du canal
title('réponse du canal')
subplot(2,1,2)
plot(f,angle(h)) % phase du canal
title('phase du canal')

% En fréquentiel, on multiplie le signal par la réponse fréquentiel du
% canal, ici on est en temporel, donc on fait le produit de convolution
xrec = conv(h,x);
figure,
subplot(2,2,1); plot(real(xrec))
title('partie réelle de x')
subplot(2,2,2); plot(imag(xrec))
title('partie imaginaire de x')
%% Ajout de bruit complexe

xrec = xrec + 0.01*(randn(size(xrec)) + j*randn(size(xrec)));
subplot(2,2,3); plot(real(xrec))
title('partie réelle de x avec bruit')
subplot(2,2,4); plot(imag(xrec))
title('partie imaginaire de x avec bruit')

%% Réception

% Ici on va tester deux structures DFVE pour recevoir le signal de départ,
% une temporelle, et une fréquentielle, la plus simple est celle temporelle
% car elle utilise des matrices triangulaire

%Matrice triangulaire contenant les coefs qui représentent les échantillons
%de la réponse du canal pendant l'emission d'un symbole OFDM, mais on les
%mets faux pour voir si notre algo LMS estime bien le canal
H0=[h(1)+0.1 0 0 0;h(2)+0.4 h(1) 0 0;h(3) h(2)+0.05 h(1) 0 ; h(4) h(3)+0.2 h(2) h(1)];
H1=[0 h(3)+0.1 h(2) h(1); 0 0 h(3)+0.5 h(2)+0.005 ; 0 0 0 h(3) ; 0 0 0 0];

%Ceci est notre matrice de TFD
w=exp(-2*pi*j/4);
MatriceTFD=[1 1 1 1;1 w w^2 w^3;1 w^2 w^4 w^6;1 w^3 w^6 w^9];
%Matrice optimales selon le critère ZF de la relation 6 de la partie 2.2 
P0=MatriceTFD*inv(H0)*inv(MatriceTFD);
P1=-MatriceTFD*inv(H0)*H1*inv(MatriceTFD);
%Matrices de la relation 8 de la partie 2.2
Q0=inv(H0);
Q1=-H1;
%% selon la figure 1
%egalisation dfve fréquentielle avec algo LMS simplifié pour estimation des
%matrices
indNew=1;
for ind = 1:NbSymNew
% décodage du symbole ind
y=fft(xrec((ind-1)*Nb+1:ind*Nb));
%si on a les valeurs pilotes, on estime PO et P1
if mod(ind,3)==0
    Xdectemp = P0*y+P1*(1./(Xdec((indNew-2)*Nb+1:(indNew-1)*Nb)))';
    %estimation des matrices
 P0=P0-0.01.*(Xdectemp-pilote)*(y)';
 P1=P1-0.01.*(Xdectemp-pilote)*(pilote)';
end
if indNew==1
    Xdec((indNew-1)*Nb+1:indNew*Nb) = P0*y;
    indNew=indNew+1;
else
    Xdec((indNew-1)*Nb+1:indNew*Nb) = P0*y+P1*(1./(Xdec((indNew-2)*Nb+1:(indNew-1)*Nb))');
    indNew=indNew+1;
end

end

%% On enleve les pilotes pour recomposer le signal
RecomposeSignal=reshape(Xdec,4,length(X)/4);
g=size(RecomposeSignal)
SignalDepartEstime=[];
for t=1:g(2)
    if(mod(t,3)==0)
    else
        SignalDepartEstime=[SignalDepartEstime;RecomposeSignal(:,t).'];
    end

end

h=size(SignalDepartEstime);
for v=1:h(2)
    for vv=1:h(1)
        SignalDepartEstime2((vv-1)*4+v)=SignalDepartEstime(vv,v);
    end
end

scatterplot(SignalDepartEstime2)
title('Constelations après structure d egalisation fréquentielle DFVE')

%démodulation
Xdemod = pskdemod(SignalDepartEstime2, M)';

erreur=1-length(find(Xdemod==R(1:356)))/length(R(1:356));




%% selon la figure 2
%Matrices de la relation 8 de la partie 2.2
for ind = 1:NbSym

if mod(ind,3)==0
    Xdectemp = fft(Q0*(xrec((ind-1)*Nb+1:ind*Nb)+Q1*(ifft((1./(Xdec2((ind-2)*Nb+1:(ind-1)*Nb))')))));
    %estimation des matrices
  Q0=Q0-0.01.*(ifft(Xdectemp)-piloteIFFT)*(xrec((ind-1)*Nb+1:ind*Nb))';
  Q1=Q1-0.01.*(ifft(Xdectemp)-piloteIFFT)*(piloteIFFT)';
end

if ind==1
    Xdec2((ind-1)*Nb+1:ind*Nb) = fft(Q0*xrec((ind-1)*Nb+1:ind*Nb));
else
    Xdec2((ind-1)*Nb+1:ind*Nb) = fft(Q0*(xrec((ind-1)*Nb+1:ind*Nb)+Q1*(ifft((1./(Xdec2((ind-2)*Nb+1:(ind-1)*Nb))'))))) ;
end
end

scatterplot(Xdec2)
title('Constelations après structure d egalisation tempo DFVE et avec pilote')

%% On enleve les pilotes pour recomposer le signal
RecomposeSignal=reshape(Xdec2,4,length(Xdec2)/4);
g=size(RecomposeSignal)
SignalDepartEstime=[];
for t=1:g(2)
    if(mod(t,3)==0)
    else
        SignalDepartEstime=[SignalDepartEstime;RecomposeSignal(:,t).'];
    end

end

h=size(SignalDepartEstime);
for v=1:h(2)
    for vv=1:h(1)
        SignalDepartEstime2((vv-1)*4+v)=SignalDepartEstime(vv,v);
    end
end

scatterplot(SignalDepartEstime2)
title('Constelations après structure d egalisation temporelle DFVE')

%démodulation
Xdemod2 = pskdemod(SignalDepartEstime2, M)';

erreur2=1-length(find(Xdemod2==R(1:356)))/length(R(1:356));
