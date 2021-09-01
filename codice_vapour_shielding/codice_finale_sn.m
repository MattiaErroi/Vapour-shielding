%Codice finale ----- parete + vapore + plasma: vapour shielding,
%evaporazione del LM, ebollizione nucleata, proprietà dipendenti dalla temperatura

clear all
close all
clc

set(0,'DefaultAxesLineWidth',0.7)
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',10)

%% ipotesi
%PARETE
%ipotesi 1): proprietà non dipendenti dalla pressione;  
%ipotesi 2): flusso imposto minore del flusso critico, quindi sempre %ebollizione nucleata
%ipotesi 3): deflusso dell'acqua completamente sviluppato e a regime (invece la temperatura della parete evolve nel tempo)
%ipotesi 4): proprietà calcolate a Twater (e non alla temperatura media di
%uscita ed entrata)
%ipotesi 5): fattore f=0.5 dal Roccella per la media pesata delle proprietà
%della CPS
%ipotesi 6): no resistenza di contatto e no resistenza tra CuCrZr e tubo (solo convettiva)
%ipotesi 7): no interferenza tra i tubi e no shape factor del 2D

%modello 1D composto da CPS e CuCrZr. Parametri geometrici e
%termofluidodinamici dal paper di Roccella. La quotatura del pezzo, non
%ricavabile interamente dal paper, è stata completata attraverso la proporzione delle misure della figura sul paper 

%VAPORE
%1) tau e Lcloud parametri liberi

%modello 0D che descrive un bilancio di particelle

%PLASMA
%1) n_e_t calcolata a partire da n_e_u, Tet, Teu ipotizzando la conservazione della pressione
%2) peso atomico della miscela calcolato come media del peso atomico di
%deuterio e trizio
%3) la funzione di radiazione non dipende da tau. I dati sono stati presi
%per tau=1ms

%modello di Lengyel

%% Parete - proprietà termofisiche
Tw=140; %[°C]
Tin=Tw;  
LCuCrZr=8.46e-3; %[m]
LCPS=2e-3; %[m] 
Ltot=LCPS+LCuCrZr;
DD=8e-3; %[m]
vel=12; %[m/s]
delta=0.8e-3; %[m]   %twisted tape thickness (fonte: SOLPS-ITER simplified heat transfer model for plasma facing components, Stefano Carli)
diameter=4*(pi*(DD^2)/4-delta*DD)/(pi*DD+2*DD-2*delta);  %diametro idraulico (fonte: SOLPS-ITER simplified heat transfer model for plasma facing components, Stefano Carli)
speed=vel*(1+pi^2/4/4)^0.5; %velocità con twisted tape (fonte: SOLPS-ITER simplified heat transfer model for plasma facing components, Stefano Carli)
RR=8.314;  %costante dei gas

%acqua (saturated water @ Tw da Incropera); b:410K, t=420K:
cp_b=4.278e3; %[J/kg/K]
cp_t=4.302e3; %[J/kg/K]
rho_b=1/1.077*1e3; %[kg/m^3]
rho_t=1/1.088*1e3; %[kg/m^3]
mu_b=200e-6; %[N*s/m^2]
mu_t=185e-6; %[N*s/m^2]
kk_b=688e-3; %[W/m/K]
kk_t=688e-3; %[W/m/K]
Pr_b=1.24;
Pr_t=1.16; 

mu_f_water=@(T) 2.4638e-5*exp(4.42e-4*50+(4.703e3-50*0.9565)/RR/(T+273.15-140.3-50*1.24e-2));  %[Pa*s] %(fonte:Dependence of Water Viscosity on Temperature and Pressure, E. R. Likhachev)

cp_w=cp_b+(Tw+273.15-410)/(420-410)*(cp_t-cp_b); %[J/kg/K]
rho_w=rho_b+(Tw+273.15-410)/(420-410)*(rho_t-rho_b); %[kg/m^3]
mu_w=mu_b+(Tw+273.15-410)/(420-410)*(mu_t-mu_b); %[N*s/m^2]
kk_w=kk_b+(Tw+273.15-410)/(420-410)*(kk_t-kk_b); %[W/m/K]
Pr_w=Pr_b+(Tw+273.15-410)/(420-410)*(Pr_t-Pr_b);

%CuCrZr-da file forniti
rho_cucrzr=@(T) 8900.*(1 - 0.000003.*(7.20e-9.*T.^3 - 9.05e-6.*T.^2+6.24e-3.*T+1.66*1e1).*(T - 20)); %T in °C %[kg/m^3]
cp_cucrzr=@(T) 6.32e-6.*T.^2 + 9.49e-2.*T + 3.88e2;  %T in °C  %[J/kg/K]
kk_cucrzr=@(T) 2.11e-7.*T.^3-2.83e-4.*T.^2 + 1.38e-1.*T + 3.23e2;  %T in °C %[W/m/K]

%Sn-da file forniti
rho_sn=@(T) 6979-0.652.*(T+273.15-505.08); %T in °C  %[kg/m^3]
mu_sn=@(T) 1e-3*10.^(-0.408+343.4./(T+273.15));  %T in °C %[N*s/m^2]
kk_sn=@(T) 13.90+0.02868.*(T+273.15);  %T in °C %[W/m/K]
cp_sn=@(T) (9.97-9.15e-3.*(T+273.15)+6.5e-6.*(T+273.15).^2)/(118.71*1e-3)*4.184;  %T in °C  %[J/kg/K]

%W-da file forniti
rho_tung=@(T) 1e3.*(19.3027-2.3786e-4.*T-2.2448e-8.*T.^2); %T in °C  %[kg/m^3]
cp_tung=@(T) 128.308+3.2797e-2.*T-3.4097e-6.*T.^2; %T in °C  %[J/kg/K]
kk_tung=@(T) 174.9274-0.1067.*T+5.0067e-5.*T.^2-7.8349e-9.*T.^3; %T in °C %[W/m/K]

%CPS-proprietà ricavate attraverso una media pesata come dal paper fornito
%(Power handling of a liquid-metal based CPS structure under high steady-state heat and particle fluxes, Morgan et al.)
ff=0.5; 
kk_CPS_sn=@(T) ff*kk_sn(T)+(1-ff)*kk_tung(T); %[W/m/K]
rho_CPS_sn=@(T) ff*rho_sn(T)+(1-ff)*rho_tung(T);  %[kg/m^3]
cp_CPS_sn=@(T) ff*cp_sn(T)+(1-ff)*cp_tung(T); %[J/kg/K]

%% Parete - calcolo dei parametri per il coefficiente di scambio termico (HTC)
Re_w=rho_w*speed*diameter/mu_w;
pressure=5; %[MPa]
T_sat=263.8; %[°C] temperatura di saturazione a 5 MPa
factor=1.15; %fattore di correzione per twisted tape (fonte: SOLPS-ITER simplified heat transfer model for plasma facing components, Stefano Carli) 
hh_st=@(Twall) factor*kk_w/diameter*0.027*(Re_w^(4/5))*(Pr_w^(1/3))*(mu_w/mu_f_water(Twall))^0.14; %correlazione di Sieder-Tate modificata %[W/m^2/K]

%% evaporazione LM
eta=1.66;
fredep=0.9;
molecular_weight=118.71*1.66054*1e-27; %[kg] 
Boltzmann=1.38*1e-23; %[J/K]
Avogadro=6.022*1e23; %[mol^-1]
enthalpy=@(T) -1285.372+28.4512*(T+273.15); %[J/mol]  %da proprietà Sn
vapor_pressure=@(T) 101325*10^(5.262-15332/(T+273.15));  %[Pa]  %da proprietà Sn
molar_flux=@(T) eta*fredep*vapor_pressure(T)/sqrt(2*pi*molecular_weight*Boltzmann*(273.15+T))/Avogadro;  %[mol/m^2/s]

%% Parete - discretizzazione
%Discretizzazione nello spazio
dx=1e-5;
xx1=(0:dx:LCuCrZr)'; %pedice1=CuCrZr
ni=length(xx1); %nodo di interfaccia
xx2=(LCuCrZr+dx:dx:Ltot)'; %pedice2=CPS
xx=[xx1;xx2];
nn=length(xx);

%Discretizzazione nel tempo
dt=1e-1;

%Preallocazione matrice dei coefficienti
AA=zeros(nn,nn);

%% Vapore - funzione analitica
tau=1e-3; %[s]
Lcloud=2e-2; %[m]

nz_fun=@(tt,Nev) Nev*tau/Lcloud*(1-exp(-tt/tau));  %[m^-3] %concentrazione iniziale del vapore nulla

%% Plasma
%dati
gamma=7;  %da Stangeby, The Plasma Boundary of Magnetic Fusion Devices
kke_par0=2390; %[W/m/eV^7/2]  %da Stangeby, The Plasma Boundary of Magnetic Fusion Devices
q_par_u=1.6e9; %[W/m^2]  %da design DTT
nn_e_u=1.3e20; %[1/m^3]  %da design DTT
LL_par=18; %[m]  %da design DTT
mm=(3.01605+2.01309)*1.66054*1e-27/2;   %peso atomico (media D e T) [kg]
beta=pi/6; %[°] %inclinazione delle linee di campo rispetto al target, da design DTT
B_theta_B_u=1.7/6.2; %(B_theta/B)_u  %da design DTT
B_theta_B_u___B_theta_B_t=3;  %(B_theta/B)_u/%(B_theta/B)_t  %da design DTT
rapp_aree=sin(beta)*B_theta_B_u/B_theta_B_u___B_theta_B_t;

%coefficiente di mitigazione per il gas nobile
coeff_mitigazione=1.8; %!!!!!! usare questo coefficiente per vedere % stato stazionario oscillatorio
% coeff_mitigazione=4.2;  %!!!!!! usare questo coefficiente per vedere % stato stazionario asintotico

%lettura funzione di radiazione
matrice_dati=load('Sn_dati.dat');  %lettura dati della funzione di radiazione per Sn
vector=reshape(matrice_dati',[],1); %sistemo in un vettore
temp=vector(1:200);  %temperature
lz_fun=vector(201:end); %valore funzione di radiazione

%definizione funzioni
nn_e_t=@(Tet,Teu) nn_e_u*Teu/Tet;  %conservazione pressione
intervallo=@(xx,yy) (linspace(xx,yy,30000))';  %definisco una funzione che genera un intervallo dati gli estremi dell'intervallo in input
rad_fun=@(xx,yy) spline(temp,lz_fun,intervallo(xx,yy)).*(intervallo(xx,yy)>=1); %interpolo i dati in base all'intervallo che viene generato. La funzione vale 0 per Te<1 eV

%iterazioni
guess=[1e9;1;1000]; %guess iniziale 
guess_2=[1e7;0.01;150]; %guess iniziale 2---> per velocizzare la risoluzione del sistema nonlineare per il plasma conviene definire un nuovo guess iniziale, più vicino alla soluzione 
toll=1e-5;
max_iter=500;

%% Condizione iniziale
figure(1)
Tm=Tin*ones(nn,1);
time(1)=0;
flux(1)=q_par_u*rapp_aree/coeff_mitigazione;   %il flusso al target si ricava da q//,u attraverso il rapporto tra le aree ed è attenuato dall'effetto del gas nobile 
Tsurf(1)=Tm(end);
Tpeak=Tw;
Tprofile=Tm;

plot(xx*1e3,Tm,'displayname',['t_0'])
hold on
grid minor
xlim([0,xx(end)*1e3])
xlabel("Spessore [mm]")
ylabel("Temperatura [°C]")
title('Andamento temperatura del target')
set(gcf,'position',[150,150,1000,570])

%% Soluzione
%Metodo Eulero Implicito con Frozen Coefficients
ii=1;
precisione=1;
precisione_1=1;
while precisione>1e-3 && precisione_1>1e-3    %per arrivare alla periodicità/stato asintotico  
    %Matrice
    alpha=(kk_cucrzr(Tprofile)./(rho_cucrzr(Tprofile).*cp_cucrzr(Tprofile))).*(xx<=LCuCrZr)+(kk_CPS_sn(Tprofile)./(rho_CPS_sn(Tprofile).*cp_CPS_sn(Tprofile))).*(xx>LCuCrZr); %CPS con Sn
    aa=alpha*dt./dx^2;
    inferiore=[-aa(2:end);0];
    main=(1+2*aa);
    superiore=[0;-aa(1:end-1)];
    DD=[inferiore,main,superiore];
    AA=spdiags(DD,-1:1,nn,nn);

    %Calcolo del HTC-metodo iterativo con Newton
    tolerance=1e-6;
    funct=@(xx) hh_st(Tm(1))*(xx-Tw)-15500*((pressure^(1.156))*(1.8*abs(xx-T_sat))^(2.046/(pressure^0.0234))); %correlazione di Bergles-Rohsenow
    dfunct=@(xx) hh_st(Tm(1))-15500*(pressure^(1.156))*(2.046/(pressure^0.0234))*((1.8*abs(xx-T_sat))^(2.046/(pressure^0.0234)))*1.8;

    [T_ONB]=myNewton(funct,dfunct,T_sat+1,tolerance); %+1 perché la funzione non è definita in T_sat  

    if Tm(1)>T_ONB   %si ha ebollizione nucleata
        qq_st=hh_st(Tm(1))*(Tm(1)-Tw);
        qq_nb=1e6*(exp(pressure/8.7)*(Tm(1)-T_sat)/22.65)^2.8;  %correlazione di Thoms-CEA 
        qq_0=1e6*(exp(pressure/8.7)*(T_ONB-T_sat)/22.65)^2.8;

        qq_tot=sqrt(qq_st^2+(qq_nb^2)*((1-qq_0/qq_nb)^2));

        hh=qq_tot/(Tm(1)-Tw);
    else
        hh=hh_st(Tm(1));    %convezione con acqua monofase sottoraffreddata
    end

    %Condizioni al contorno

        %Robin primo nodo
        AA(1,1)=-(1+dx*hh/kk_cucrzr(Tprofile(1)));
        AA(1,2)=1;

        %Interfaccia
        AA(ni,ni-1)=kk_cucrzr(Tprofile(ni-1))/kk_CPS_sn(Tprofile(ni-1));  
        AA(ni,ni)=-1-kk_cucrzr(Tprofile(ni))/kk_CPS_sn(Tprofile(ni));
        AA(ni,ni+1)=1;

        %Neumann non omogenea
        AA(end,end-1)=-1;
        AA(end,end)=1;

        %Vettore termini noti
        bb=Tm;
        bb(1)=-Tw*hh*dx/kk_cucrzr(Tprofile(1));
        bb(ni)=0;
        bb(end)=(flux(ii)-molar_flux(Tm(end))*enthalpy(Tm(end)))*dx/kk_CPS_sn(Tprofile(end));        
        
        %Soluzione
        TT=AA\bb;
        
        if max(TT)>Tpeak
            Tpeak=max(TT);
            flag=1;
        else
            flag=0;
        end
                        
        if flag==0 
         precisione=norm(max(TT)-Tpeak)/norm(Tpeak-Tw);   %differenza tra i massimi di ogni periodo
        end
        
        precisione_1=norm(TT-Tm)/norm(TT-Tw);  %per lo stato asintotico
        
        ii=ii+1;
        time(ii)=time(ii-1)+dt;
        
        Tsurf(ii)=TT(end);
        Tprofile=Tm;   %Tprofile è il profilo di temperatura all'istante precedente, ma che a differenza di Tm non viene aggiornato. È preferibile valutare le proprietà a Tprofile e non a Tm per i frozen coefficients perché si ha periodicità
        
        %Calcolo concentrazione vapore
        nz=nz_fun(time(ii),molar_flux(TT(end))*Avogadro);
        
        %Plasma-soluzione del sistema nonlineare con il metodo di Newton
        qq_par=@(qq,Tet,Teu,TT) (abs(qq^2+2*kke_par0*nz/nn_e_t(Tet,Teu)*nn_e_u^2*Teu^2*trapz(TT,sqrt(TT).*spline(temp,lz_fun,TT).*(TT>=1))))^(1/2);  %funzione ausiliaria

        funx1=@(qq,Tet,Teu) qq-(abs(q_par_u^2-2*kke_par0*nz/nn_e_t(Tet,Teu)*nn_e_u^2*Teu^2*trapz(intervallo(Tet,Teu),sqrt(intervallo(Tet,Teu)).*rad_fun(Tet,Teu))))^(1/2);
        funx2=@(qq,Tet,Teu) qq-gamma*nn_e_u*Teu*(1.602e-19)/2*sqrt(2*Tet*(1.602e-19)/mm);
        funx3=@(qq,Tet,Teu) LL_par-kke_par0*trapz(intervallo(Tet,Teu),intervallo(Tet,Teu).^(5/2)./qq_par(qq,Tet,Teu,intervallo(Tet,Teu)));

        fun2a=@(x) [funx1(x(1),x(2),x(3));funx2(x(1),x(2),x(3));funx3(x(1),x(2),x(3))];  %vettore avente per righe le equazioni del sistema
        jfun2a=@(x) numerical_jacobian(fun2a,x,1e-6);
        
        if nz>1e21  %conviene usare guess_2
          [yy,err,residual,niter]=myNewton_Jac(fun2a,jfun2a,guess_2,toll,max_iter);
        else
          [yy,err,residual,niter]=myNewton_Jac(fun2a,jfun2a,guess,toll,max_iter); 
        end
        
        %Aggiornamento
        flux(ii)=yy(1)*rapp_aree/coeff_mitigazione;   %passo da q//,t a q''t
        Tm=TT;
        
        %Plot
        plot(xx*1e3,TT,'displayname',['t=' num2str(time(ii)) 's']);
        pause(0.01)
        legend('location','westoutside')
end

figure(2)
plot(time,flux*1e-6,'r')
grid minor
xlim([0,time(end)])
xlabel("Tempo [s]")
ylabel("Flusso termico [MW/m^2]")
title('Flusso al target q''''_t nel tempo')

figure(3)
plot(time,Tsurf)
grid minor
xlabel("Tempo [s]")
ylabel("Temperatura [K]")
xlim([0,time(end)])
title('Temperatura all''estremo alto del dominio nel tempo')

figure(4)
plot(xx*1e3,Tin*ones(size(xx)),'displayname',['t_0'])
hold on
plot(xx*1e3,Tprofile,'displayname',['t=' num2str(time(end-1)) 's'])
plot(xx*1e3,TT,'displayname',['t=' num2str(time(end)) 's'])
hold on
grid minor
xlim([0,xx(end)*1e3])
xlabel("Spessore [mm]")
ylabel("Temperatura [°C]")
title('Profili di temperatura del target')
legend('location','best')

%% andamento del HTC
temp=(130:310);
for ii=1:length(temp)
    funct=@(xx) hh_st(temp(ii))*(xx-Tw)-15500*((pressure^(1.156))*(1.8*abs(xx-T_sat))^(2.046/(pressure^0.0234)));
    dfunct=@(xx) hh_st(temp(ii))-15500*(pressure^(1.156))*(2.046/(pressure^0.0234))*((1.8*abs(xx-T_sat))^(2.046/(pressure^0.0234)))*1.8;

    [T_ONB]=myNewton(funct,dfunct,T_sat+1,toll); 
    
    if temp(ii)>T_ONB
        qq_st=hh_st(temp(ii))*(temp(ii)-Tw);
        qq_nb=1e6*(exp(pressure/8.7)*(temp(ii)-T_sat)/22.65)^2.8;    %Thoms-CEA correlation
        qq_0=1e6*(exp(pressure/8.7)*(T_ONB-T_sat)/22.65)^2.8;

        qq_tot=sqrt(qq_st^2+(qq_nb^2)*((1-qq_0/qq_nb)^2));

        hh=qq_tot/(temp(ii)-Tw);
    else
        hh=hh_st(temp(ii));
    end
    HTC(ii)=hh;
end
    
figure(5)
plot(temp,HTC*1e-3)
grid minor
xlim([temp(1),temp(end)])
xlabel("Temperatura [°C]")
ylabel("HTC [kW/m^2/K]")
title('HTC in funzione della temperatura interna del tubo')

%% Calcolo del flusso critico - correlazione di Tong75 modificata per twisted tape (da design ITER, fonte: SOLPS-ITER simplified heat transfer model for plasma facing components, Stefano Carli)
f_0=8/(Re_w^0.6)*(diameter/(12.7*1e-3))^0.32;
i_fg=1639.39e3; %calore latente di vaporizzazione a 5MPa [J/kg]
rho_v=26.67; %densità vapore a Tsat [kg/m^3]
J_a=cp_w*(T_sat-Tw)/i_fg*rho_w/rho_v;
Cf=1.67; %fattore di correzione - da design ITER, fonte: SOLPS-ITER simplified heat transfer model for plasma facing components, Stefano Carli

critical_flux=Cf*0.23*f_0*rho_w*speed*i_fg*(1+0.00216*((pressure/22.09)^1.8)*J_a*Re_w^0.5);

fprintf('Il flusso critico è %.2f [MW/m^2]\n',critical_flux*1e-6)
fprintf('Il flusso iniziale al target è %.2f [MW/m^2]\n',flux(1)*1e-6)


