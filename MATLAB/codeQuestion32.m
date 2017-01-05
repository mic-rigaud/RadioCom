close all, clear all, clc;
%% Création du signal OFDM

% Nombre d'états de la QAM.
M = 4;
% Nombre de porteuses dans le symbole OFDM
Nb = 4;
%Nombre de symboles OFDM dans la simulation par fréquence
NbSym = 10;
%Nombre total de symbole
NbSymTot=Nb*NbSym;
% Tirage aléatoire d'entiers allant de 0 à M-1, càd un symbole pour un état
R = randint(NbSymTot,1,M);
% Mise en constellation QAM. Quand 4 états, 4 constellation : 1 -1 j et -j
X = pskmod(R, M) ;scatterplot(X);title('Constelations du signal à l emission après modulation pi/4-QPSK');
 figure
% Création signal OFDM
x = zeros(size(X));
for ind = 1:NbSym
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
%de la réponse du canal pendant l'emission d'un symbole OFDM
H0=[h(1) 0 0 0;h(2) h(1) 0 0;h(3) h(2) h(1) 0 ; h(4) h(3) h(2) h(1)];
H1=[0 h(3) h(2) h(1); 0 0 h(3) h(2) ; 0 0 0 h(3) ; 0 0 0 0];
%Ceci est notre matrice de TFD
w=exp(-2*pi*j/4);
MatriceTFD=[1 1 1 1;1 w w^2 w^3;1 w^2 w^4 w^6;1 w^3 w^6 w^9];

%Matrice optimales selon le critère ZF de la relation 6 de la partie 2.2 
P0=MatriceTFD*inv(H0)*inv(MatriceTFD);
P1=-MatriceTFD*inv(H0)*H1*inv(MatriceTFD);

%% selon la figure 1
for ind = 1:NbSym
% décodage du symbole ind
y=fft(xrec((ind-1)*Nb+1:ind*Nb));
% sauvegarde du ième symbole décodé
if ind==1
    Xdec((ind-1)*Nb+1:ind*Nb) = P0*y;
else
    Xdec((ind-1)*Nb+1:ind*Nb) = P0*y+P1*(1./(Xdec((ind-2)*Nb+1:(ind-1)*Nb)))';
end
end
scatterplot(Xdec)
title('Constelations après structure d egalisation fréquentielle DFVE')

%démodulation
X = pskdemod(Xdec, M)';

%calcul des erreurs
erreur=1-length(find(X==R))/length(R);

%% selon la figure 2
%Matrices de la relation 8 de la partie 2.2
Q0=inv(H0);
Q1=-H1;
for ind = 1:NbSym
% décodage du symbole ind
y=fft(xrec((ind-1)*Nb+1:ind*Nb));
% sauvegarde du ième symbole décodé
if ind==1
    Xdec2((ind-1)*Nb+1:ind*Nb) = fft(Q0*xrec((ind-1)*Nb+1:ind*Nb));
else
    Xdec2((ind-1)*Nb+1:ind*Nb) = fft(Q0*(xrec((ind-1)*Nb+1:ind*Nb)+Q1*(ifft((1./(Xdec2((ind-2)*Nb+1:(ind-1)*Nb))'))))) ;
end
end
scatterplot(Xdec2)
title('Constelations après structure d egalisation temporelle DFVE')
%démodulation
X2 = pskdemod(Xdec2, M)';
erreur2=1-length(find(X2==R))/length(R);