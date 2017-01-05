close all, clear all, clc;
%% Cr�ation du signal OFDM

% Nombre d'�tats de la QAM.
M = 4;
% Nombre de porteuses dans le symbole OFDM
Nb = 4;
%Nombre de symboles OFDM dans la simulation par fr�quence
NbSym = 10;
%Nombre total de symbole
NbSymTot=Nb*NbSym;
% Tirage al�atoire d'entiers allant de 0 � M-1, c�d un symbole pour un �tat
R = randint(NbSymTot,1,M);
% Mise en constellation QAM. Quand 4 �tats, 4 constellation : 1 -1 j et -j
X = pskmod(R, M) ;scatterplot(X);title('Constelations du signal � l emission apr�s modulation pi/4-QPSK');
 figure
% Cr�ation signal OFDM
x = zeros(size(X));
for ind = 1:NbSym
% calcul i�me symbole OFDM
symbole=ifft(X((ind-1)*Nb+1:ind*Nb));
% sauvegarde du symbole ind dans x
x((ind-1)*Nb+1:ind*Nb) = symbole; % bien comprendre que l'info x(1) est sur la fr�quence de la sous-porteuse 1, x(2) sur la 2,..., x(5) sur la 1,...
end
subplot(2,1,1); plot(real(x))
title('partie r�elle de x')
subplot(2,1,2); plot(imag(x))
title('partie imaginaire de x')



%% Canal de propagation
figure;
% passage dans le canal multi-trajet
% On a mit des fr�quences comme pour le canal 1 du wifi
% Norme 802.11a : espacement entre sous-porteuses 0.3125MHz
fcentre=2.412e9;
espacementSP=0.3125e6;
f=[fcentre,fcentre+espacementSP,fcentre+2*espacementSP,fcentre+3*espacementSP]; 
h=ones(1,4)+(0.4+j*0.2).*(1./f); % r�ponse fr�quenciel d'un canal �cho
subplot(2,1,1)
plot(f,abs(h)) % r�ponse du canal
title('r�ponse du canal')
subplot(2,1,2)
plot(f,angle(h)) % phase du canal
title('phase du canal')

% En fr�quentiel, on multiplie le signal par la r�ponse fr�quentiel du
% canal, ici on est en temporel, donc on fait le produit de convolution
xrec = conv(h,x);
figure,
subplot(2,2,1); plot(real(xrec))
title('partie r�elle de x')
subplot(2,2,2); plot(imag(xrec))
title('partie imaginaire de x')
%% Ajout de bruit complexe

xrec = xrec + 0.01*(randn(size(xrec)) + j*randn(size(xrec)));
subplot(2,2,3); plot(real(xrec))
title('partie r�elle de x avec bruit')
subplot(2,2,4); plot(imag(xrec))
title('partie imaginaire de x avec bruit')
%% R�ception

% Ici on va tester deux structures DFVE pour recevoir le signal de d�part,
% une temporelle, et une fr�quentielle, la plus simple est celle temporelle
% car elle utilise des matrices triangulaire

%Matrice triangulaire contenant les coefs qui repr�sentent les �chantillons
%de la r�ponse du canal pendant l'emission d'un symbole OFDM
H0=[h(1) 0 0 0;h(2) h(1) 0 0;h(3) h(2) h(1) 0 ; h(4) h(3) h(2) h(1)];
H1=[0 h(3) h(2) h(1); 0 0 h(3) h(2) ; 0 0 0 h(3) ; 0 0 0 0];
%Ceci est notre matrice de TFD
w=exp(-2*pi*j/4);
MatriceTFD=[1 1 1 1;1 w w^2 w^3;1 w^2 w^4 w^6;1 w^3 w^6 w^9];

%Matrice optimales selon le crit�re ZF de la relation 6 de la partie 2.2 
P0=MatriceTFD*inv(H0)*inv(MatriceTFD);
P1=-MatriceTFD*inv(H0)*H1*inv(MatriceTFD);

%% selon la figure 1
for ind = 1:NbSym
% d�codage du symbole ind
y=fft(xrec((ind-1)*Nb+1:ind*Nb));
% sauvegarde du i�me symbole d�cod�
if ind==1
    Xdec((ind-1)*Nb+1:ind*Nb) = P0*y;
else
    Xdec((ind-1)*Nb+1:ind*Nb) = P0*y+P1*(1./(Xdec((ind-2)*Nb+1:(ind-1)*Nb)))';
end
end
scatterplot(Xdec)
title('Constelations apr�s structure d egalisation fr�quentielle DFVE')

%d�modulation
X = pskdemod(Xdec, M)';

%calcul des erreurs
erreur=1-length(find(X==R))/length(R);

%% selon la figure 2
%Matrices de la relation 8 de la partie 2.2
Q0=inv(H0);
Q1=-H1;
for ind = 1:NbSym
% d�codage du symbole ind
y=fft(xrec((ind-1)*Nb+1:ind*Nb));
% sauvegarde du i�me symbole d�cod�
if ind==1
    Xdec2((ind-1)*Nb+1:ind*Nb) = fft(Q0*xrec((ind-1)*Nb+1:ind*Nb));
else
    Xdec2((ind-1)*Nb+1:ind*Nb) = fft(Q0*(xrec((ind-1)*Nb+1:ind*Nb)+Q1*(ifft((1./(Xdec2((ind-2)*Nb+1:(ind-1)*Nb))'))))) ;
end
end
scatterplot(Xdec2)
title('Constelations apr�s structure d egalisation temporelle DFVE')
%d�modulation
X2 = pskdemod(Xdec2, M)';
erreur2=1-length(find(X2==R))/length(R);