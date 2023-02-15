clear all
close all
clc

%%%%%%%%%%%% numbering convention %%%%%%%%%%%%
    
                %%%%%   %%%%%
                % 1 %   % 2 %
                %%%%%   %%%%%
            
                %%%%%   %%%%%
                % 4 %   % 3 %
                %%%%%   %%%%%


%% Synapse Parameters

thetaE = -50; % step 0
thetaD = -50; % step 0
thetaF = -50; % step 0
thetaS = 17.5; % step 2

nu = 10; % step 0

alpha = 0.754312006335462; % step 2
beta = 0.039069399370546; % step 2

gS = 0.120679264063933; % step 4
gE = 0.005179474679231; % step 4
gDgSratio = 0.011758053907164; % step 0
gF41gF14ratio = 2.332706766917304; % step 0
gF41 = 0.004756160632129; % step 3
gF14 = gF41*gF41gF14ratio; % step 3

EsynS = -80;
EsynE = 60;

%%% burst period of single neuron (one value for each value of Ic_) used to calculate delays of D synapses %%%
period = [94.3109051428739	95.6826953384914	96.8043762328233	98.4040549166303	99.8780657331680...
        101.280159359833	102.838773224492	104.326305214116	105.973457733726	107.546295603950...
        109.110160345552	110.941767521370	112.644072120751	114.300251095713	116.329146262383...
        118.183043638007	120.057919413082	122.170807444693	124.188943347961	126.259416526871...
        128.378117009658	130.520734719829	132.967085873819	135.298145973170	137.697957134677...
        140.175527976435	142.730825225410	145.367998767864	148.097782334000	150.921366168908...
        153.846602462613	156.862571372874	160.369953495542	163.733589505284	167.260927851001...
        170.976753540141	174.904678990471	179.068338090176	183.506458007384	188.255911651234...
        193.373151664137	198.921092874419	204.979588201391	211.668548387157	219.149189289523...
        227.663957669443	237.588805589099	249.582608632000	265.202497841058	288.572356306551];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_ = period/2; % step 0

%% Neuron Parameters

gna = 100;
ENa = 50;
gk = 10;
Ek = -95;
gca = 1.75;
gl = 0.05;
El = -78;
C = 1;
Ca0 = 2;
d = 1;
KT = 0.0001;
Kd = 0.0001;
gd = 0.0001;
D = 0;
tau = 0.33;
tTm1 = 1;
tTm2 = 1;
tTh1 = 1;
tTh2 = 1;

N=4; % number of neurons in the network

%% Gait selection

Ic_=linspace(-0.43,0.13,50);
gait = 'bound';
switch gait
    case 'bound'
        Ic12 = Ic_(3);
        Ic34 = Ic_(8); 
        delta12 = delta_(3);
        delta34 = delta_(8); 
    case 'trot'
        Ic12 = Ic_(37);
        Ic34 = Ic_(39); 
        delta12 = delta_(37);
        delta34 = delta_(39); 
    case 'walk'
        Ic12 = Ic_(49);
        Ic34 = Ic_(44);
        delta12 = delta_(49);
        delta34 = delta_(44);
end

%% Network building

% creating synapses objects
dummy = FTM_synapse_model(); % empty synapse corresponding to null weight, only used as place holder
E_synapse = FTM_synapse_model(nu,thetaE);
F_synapse = FTM_synapse_model(nu,thetaF);
S_synapse12 = alpha_delayed_synapse_model(alpha,beta,nu,thetaS,delta12,gDgSratio);
S_synapse34 = alpha_delayed_synapse_model(alpha,beta,nu,thetaS,delta34,gDgSratio);

% filling matrices of excitatory and inhibitory synapses
synapses_ex = {dummy, E_synapse, dummy, dummy;
                E_synapse, dummy, dummy, dummy;
                dummy, dummy, dummy, E_synapse;
                dummy, dummy, E_synapse, dummy};

synapses_in = {dummy, S_synapse12, dummy, F_synapse;
                S_synapse12, dummy, F_synapse, dummy;
                dummy, F_synapse, dummy, S_synapse34;
                F_synapse, dummy, S_synapse34, dummy};

% filling matrices of synaptic weights
g_in =  [0, gS, 0, gF14; 
         gS, 0, gF14, 0;
         0, gF41, 0, gS;
         gF41, 0, gS, 0];
     
g_ex = gE*[0, 1, 0, 0; 
           1, 0, 0, 0;
           0, 0, 0, 1;
           0, 0, 1, 0];

g_el =  zeros(4);

% creating neuron objects
neuron1 = RE_model(gna,ENa,gk,Ek,gca,gl,El,C,Ic12,Ca0,d,KT,Kd,gd,D,EsynE,tau,tTm1,tTm2,tTh1,tTh2);
neuron2 = RE_model(gna,ENa,gk,Ek,gca,gl,El,C,Ic12,Ca0,d,KT,Kd,gd,D,EsynE,tau,tTm1,tTm2,tTh1,tTh2);
neuron3 = RE_model(gna,ENa,gk,Ek,gca,gl,El,C,Ic34,Ca0,d,KT,Kd,gd,D,EsynE,tau,tTm1,tTm2,tTh1,tTh2);
neuron4 = RE_model(gna,ENa,gk,Ek,gca,gl,El,C,Ic34,Ca0,d,KT,Kd,gd,D,EsynE,tau,tTm1,tTm2,tTh1,tTh2);

% filling matrix of neurons
neurons = {neuron1,neuron2,neuron3,neuron4};

% creating network object
netw = CPG(N,neurons,g_in,g_ex,g_el,EsynS,EsynE,synapses_in,synapses_ex);

%% Simulation

% impose phase lag initial conditions, both favorable (phase lag of the selected gait)
% and limit case (close to adjacent gait/gaits)
switch gait
    case 'bound'
        phi0 = [0, 0.5, 0.5;
                0.45, 0.05, 0.45]; 
        Nphi0 = 2;
    case 'trot'
        phi0 = [0.5, 0, 0.5;
                0.05, 0.45, 0.45;
                0.45, 0.7, 0.3]; 
        Nphi0 = 3;
    case 'walk'
        phi0 = [0.5, 0.75, 0.25;
                0.45, 0.05, 0.45]; 
        Nphi0 = 2;
end

% estimate initial conditions of the networks
options1.Vth = -30;
options1.x0 = [0,0,0,0,0,0,0.000001];
options1.Ttrans = 2000;
IC = netw.getCIfromPhi(phi0,options1);

% simulate networks
options2.Vth = -30;
options2.integrator = 'odeint';
options2.integratorOptions = odeset('RelTol',1e-9,'AbsTol',1e-12);
Tspan = [0,6000]; % simulate for 6 seconds

for i=1:Nphi0
[TT1{i},XX1{i}] = netw.sim(Tspan,IC(i,:),options2);
end

%% Plot results

figure
for i=1:Nphi0
subplot(Nphi0,1,i)
hold on
plot(TT1{i},XX1{i}(:,1),'g')
plot(TT1{i},XX1{i}(:,8),'b')
plot(TT1{i},XX1{i}(:,15),'r')
plot(TT1{i},XX1{i}(:,22),'y')
legend('1','2','3','4')
xlim([Tspan(2)-1000 Tspan(2)]) % show only 1 second of steady state behavior
ylim([-80 60])
ylabel('V_i [mV]')
xlabel('t [s]')
title(gait)
end