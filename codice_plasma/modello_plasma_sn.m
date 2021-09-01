%Modello plasma di Lengyel con il metodo di Newton 

clear all
close all
clc

set(0,'DefaultAxesLineWidth',0.7)
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',10)

%ipotesi:
%1) n_e_t calcolata a partire da n_e_u, Tet, Teu ipotizzando la conservazione della pressione
%2) peso atomico della miscela calcolato come media del peso atomico di
%deuterio e trizio
%3) la funzione di radiazione è considerata indipendente da tau. I dati sono stati presi
%per tau=1ms

%% dati
gamma=7;  %da Stangeby, The Plasma Boundary of Magnetic Fusion Devices
kke_par0=2390; %[W/m/eV^7/2]  %da Stangeby, The Plasma Boundary of Magnetic Fusion Devices
q_par_u=1.6e9; %[W/m^2]  %da design DTT
nn_e_u=1.3e20; %[1/m^3]  %da design DTT
LL_par=18; %[m]  %da design DTT

nz=1e20;  %ipotesi (nel codice finale questo termine viene dal modello del vapore) [m^-3]

mm=(3.01605+2.01309)*1.66054*1e-27/2;   %peso atomico (media D e T) [kg]

%% interpolazione funzione radiazione
matrice_dati=load('Sn_dati.dat');  %lettura dati della funzione di radiazione per Sn
vector=reshape(matrice_dati',[],1); %sistemo in un vettore
temp=vector(1:200);  %temperature
lz_fun=vector(201:end); %valore funzione di radiazione

zz=linspace(temp(1),temp(end),30000);
ss=spline(temp,lz_fun,zz);

figure(1)
loglog(temp,lz_fun,'g',zz,ss,'r--')
grid on
xlabel("Te [eV]")
ylabel("P_{rad}/(n_e n_z) [Wm^3]")
title('Funzione di radiazione per Sn con tau=1ms')
legend('funzione','interpolazione','location','best')

%% soluzione-metodo di Newton (sistema nonlineare con 3 equazioni in 3 incognite)
nn_e_t=@(Tet,Teu) nn_e_u*Teu/Tet;  %conservazione pressione
intervallo=@(xx,yy) (linspace(xx,yy,30000))';  %definisco una funzione che genera un intervallo dati gli estremi dell'intervallo in input
rad_fun=@(xx,yy) spline(temp,lz_fun,intervallo(xx,yy)).*(intervallo(xx,yy)>=1); %interpolo i dati in base all'intervallo che viene generato. Per valori di Te<1eV la funzione viene posta a 0

qq_par=@(qq,Tet,Teu,TT) (abs(qq^2+2*kke_par0*nz/(nn_e_t(Tet,Teu))*nn_e_u^2*Teu^2*trapz(TT,sqrt(TT).*spline(temp,lz_fun,TT).*(TT>=1))))^(1/2);  %funzione ausiliaria

%definizione funzione

funx1=@(qq,Tet,Teu) qq-(abs(q_par_u^2-2*kke_par0*nz/(nn_e_t(Tet,Teu))*nn_e_u^2*Teu^2*trapz(intervallo(Tet,Teu),sqrt(intervallo(Tet,Teu)).*rad_fun(Tet,Teu))))^(1/2);
funx2=@(qq,Tet,Teu) qq-gamma*nn_e_u*Teu*(1.602e-19)/2*sqrt(2*Tet*(1.602e-19)/mm);
funx3=@(qq,Tet,Teu) LL_par-kke_par0*trapz(intervallo(Tet,Teu),intervallo(Tet,Teu).^(5/2)./qq_par(qq,Tet,Teu,intervallo(Tet,Teu)));

fun2a=@(x) [funx1(x(1),x(2),x(3));funx2(x(1),x(2),x(3));funx3(x(1),x(2),x(3))];  %vettore avente per righe le equazioni del sistema nonlineare da risolvere

%calcolo jacobiana numerica

jfun2a=@(x) numerical_jacobian(fun2a,x,1e-5);

%iterazioni
guess=[1e9;1;1000];  %guess iniziale
toll=1e-5;
max_iter=500;

[yy,err,residual,niter]=myNewton_Jac(fun2a,jfun2a,guess,toll,max_iter);

fprintf('q//,t=%.2f MW/m^2\n',yy(1)*1e-6)
fprintf('Te,t=%.2f eV\n',yy(2))
fprintf('Te,u=%.2f eV\n',yy(3))

