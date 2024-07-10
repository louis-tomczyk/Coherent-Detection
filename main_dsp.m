% ---------------------------------------------
% ----- INFORMATIONS -----
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-01-13
%
% ----- Main idea -----
%   Simulating TX, PROPAGATION and RX side of a complete optical
%   telecommunication system
%
%   - transmitter side
%   - propagation
%   
% ----- BIBLIOGRAPHY -----
%   Functions           : Optilux 2009 and 2022 versions
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2009 and 2022
%   Title of program    : Optilux
%   Code version        : 2009 and 2022
%   Type                : Optical simulator toolbox - source code
%   Web Address         : v2009 --- Partage Zimbra "Phd Louis Tomczyk"/
%                               Optilux/Optilux v2009
%                         v2022 --- https://optilux.sourceforge.io/
% -----------------------
%   Articles
%   Author              : Alix MAY
%   Title               : Receiver-Based Experimental Estimation of Power
%                         Losses in Optical Networks
%   Jounal              : IEEE - Photonics Technology Letters
%   Volume - N°         : 33-22
%   Date                : 2021-11-22
%   DOI                 : 10.1109/LPT.2021.3115627
% ---------------------------------------------

%% MAINTENANCE
rst
init_step()

%% PARAMETERS
        
%%% Modulation parameters
tx.modfor       = 'qpsk';   % []    {qpsk,16qam}
tx.Nch          = 1;        % []    {integer>0}
tx.PdBm         = 0;        % [dBm] power at fibre input, sum of all channel power
% tx.seed         = 1;

%%% Axis parameters
Axis.Nsymb      = 2^14;     % []    {2^(integer>0)}
Axis.symbrate   = 32;       % [Gbd] {integer>0}

%%% Laser parameters
las.n_polar     = 2;        % []    {1,2}
las.type        = 'ideal';  %       {ideal,LRL,HPFL}

%%% Optical amplifiers
% ft.type         = 'B2B';
amp.type        = 'noAWGN'; %     {noAWGN,ideal,classic}
amp.Nspan       = 1;        % []    number of spans


%%% Transmission fiber
tx              = set_tx(tx);
Axis            = set_axis(tx,Axis);
las             = set_las(tx,las);
ft              = set_ft(las,ft);
ft.length       = 100e3;
amp             = set_topology(tx,ft,amp);
dsp             = set_dsp();
rx              = set_rx(Axis,dsp,tx,las,ft,amp);

rx.CMA.method   = "phenix";      % {simu,phenix}
% {phenix -> step = 2,simu, x15 quicker, but lower success rate}
%

rx.CMA.taps     = 11;
rx.CMA.Ntrain   = Axis.Nsamp;
rx.CMA.lr_train = 1e-3;
rx.CMA.lr_eq    = 1e-3;

rx.CMA.singularity.loops    = 1;
rx.CMA.plot_final           = 0;
rx.CMA.singularity.plot     = 0;
% rx.CMA.flag.thresh          = 1/50;

inigstate(Axis.Nsamp,Axis.fs)

%% SIMULATIONS

% Fibre input field
% Link output field
%   {SSFM,RP1,CD-MIMO,only-MIMO (XX-MIMO: add Axis as input and ft as output)}
% Optical filtering, wavelegnth selection to come
% downsampling, quantization to come
% DSP filtering
% Chromatic dispersion compensation
% Polarisation demultiplexing


for k = 1:10
    Ein         = TX(Axis,las,tx);
    [Elink,ft]  = channel_prop(Ein,ft,amp,'CD-MIMO',Axis);
    EfiltOpt    = filt_opt(Elink,Axis,rx);
    [Eadc,Axis] = adc(Elink,Axis,rx);
    EfiltElec   = filt_elec(Eadc,Axis);
    Ecdc        = cdc(EfiltElec,Axis,ft,rx);
    [Emimo,rx]  = mimo(Ecdc,Axis,rx);
    Ecpc        = cpc(Emimo,tx,rx);
end

main_dsp_results

