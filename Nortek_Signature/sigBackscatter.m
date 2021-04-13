function Data = sigBackscatter( Data, Config , mode )
% SIGBACKSCATTER estimates the backscatter coefficient from beam amplitudes
% for Signature-series ADCPs
%
%   Data = sigBackscatter( Data, Config ) calculates beam backscatter and
%   adds Sv variables for the Average mode data in the structure 'Data'.
%
%   Data = sigBackscatter( Data, Config , mode ) allows specification of
%   the input data mode as 'avg' or 'burst' (corresponding to Average or
%   Burst structure variables).  The function can act on multiple data
%   types by including different modes by including a cell array of modes:
%   e.g. {'avg','burst'}
%
%
%   Notes:  
%   (1) This function is developed to operate on Data structures that are
%   output by converting raw .ad2cp data to .mat files using MIDAS
%   software.  Data converted with Signature Deployment software may not
%   have matching variable names.
%   (2) From previous experience, the recorded temperature and calculated
%   speed-of-sound given by the instruments is slightly wrong.  I have run
%   other functions to correct these issues and simultaneously re-create a
%   salinity data vector.  This funciton checks for those corrected values
%   (using the variable names I've assigned), and if they don't exist it
%   will use default values from the Data and Config structures.  
%
%   S.D.Brenner, 2019

%% Parse inputs

    if nargin < 3 || isempty(mode); mode = 'avg'; end

    
    % Parse mode choice
    %   ( Note, 'mode' options could have instead been the 'dataWordChoices'
    %     values, but instead are 'modeChoices' to be consistent with other
    %     Nortek and Signature codes)
    modeChoices = {'avg','burst'};
    dataWordChoices = {'Average','AverageIce','Burst'};
    [modeLog,modeInd] = ismember( lower(mode) , modeChoices );
    if ~modeLog
        error('The input variable ''mode'' must be one of: ''avg'', ''ice'', or ''burst''');
    elseif length(modeLog)>1
        % If multiple mode words are entered, recursively run this script for
        % each of the individually (this may break something)
        for n = 1:length(modeLog)
            modeN = modeChoices{modeInd(n)};
            Data = sigEnsembleAvg(Data,modeN,flds,func);
        end
        return;
    else
        dataModeWord = dataWordChoices{modeInd};
    end

%% Extract data from structure

f = Config.plan_frequency;              % [kHz] ADCP tranmission frequency
D = Data.([dataModeWord,'_Pressure']);  % [dbar] (~[m] depth)
theta = Config.beamConfiguration1_theta;
r = Data.([dataModeWord,'_Range']) / cosd(theta);
for n = 1:4
    E(:,:,n) = Data.([dataModeWord,'_AmpBeam',num2str(n)]);
end

% check for "corrected" sound speed
if isfield( Data, [dataModeWord,'_SpeedOfSound_Corrected'] )
    cs = Data.([dataModeWord,'_SpeedOfSound_Corrected']);
else
    cs = Data.([dataModeWord,'_SpeedOfSound']);
end


% check for "corrected" water temperature
if isfield( Data, [dataModeWord,'_WaterTemperature_Corrected'] )
    T = Data.([dataModeWord,'_WaterTemperature_Corrected']);
else
    T = Data.([dataModeWord,'_WaterTemperature']);
end
TK = T + 273.15; % [K] temperature in kelvin


% check for re-constructed salinity
if isfield( Data, [dataModeWord,'_Salinity_Reconstructed'] )
    S = Data.([dataModeWord,'_Salinity_Reconstructed']);
else
    S = Config.salinity;
end



% check for re-constructed pH (I don't have this)
if isfield( Data, [dataModeWord,'_pH_Reconstructed'] )
    pH = Data.([dataModeWord,'_pH_Reconstructed']);
else
    pH = 7.5;
end


%% Calculate transmission loss
% Transmission loss depends on the absorption coefficient, 'alpha'.  The
% following uses the Francois and Garrison, 1982 method for estimatin alpha
% (see doi: 10.1121/1.388673)

% Boric acid contribution
A1 = (8.86./cs) * 10.^( 0.78*pH-5 );
P1 = 1;
f1 = 2.8 * (S/35).^0.5 .* 10.^( 4 - 1245./TK );
alfBA = A1 .* P1 .* f1 .* f.^2 ./( f.^2 + f1.^2);

% Magnesium sulfide contribution
A2 = 21.44*(S./cs).*( 1 + 0.025*T );
P2 = 1 - 1.37e-4*D + 6.2e-9*D.^2;
f2 = (8.17 * 10.^(8-1990./TK) )./( 1 + 0.0018*(S-35) );
alfMS = A2 .* P2 .* f2 .* f.^2 ./( f.^2 + f2.^2);

% Pure water contribution
A3 = 4.973e-4 - 2.59e-5 * T +  9.11e-7 * T.^2 - 1.50e-8 * T.^3;
P3 = 1 - 3.83e-5 * D + 4.9e-10 * D.^2; 
alfPW = A3 .* P3 .* f.^2;

% Total sound absorption
alpha = alfBA + alfMS + alfPW ; % [ dB/km ]
alpha = alpha/1000;               % convert from [ dB/km ] to [ dB/m ]

% Calculate 2-way transmission loss:
[R,Alf] = meshgrid(r,alpha);
TL = 20*log10(R) + 2*Alf.*R;

%% Calculate relative backscatter 

% Calculate transducer echo levels (RL)
% ( From Dwinovantyo et. al., 2017; doi.org/10.1155/2017/4890421 )
Te = Data.([dataModeWord,'_MagnetometerTemperature']) + 273.15; % [K] temperature of electronics
Kc = 127.3./Te;     % Dwinovantyo et. al., 2017, eqn 6
Kc = 1;
[N,M,L] = size(E);
Kc = repmat(Kc,[1,M,L]);
Er = min(E(:)); % [dB] instrument noise (found from histogram of all echo data)
% RL = Kc.*(E-Er);
RL = 10*log10( 10.^((Kc.*E)/10) - 10.^((Kc.*Er)/10) );  % Gostiuax2010 eqn 5
    
% Calculate relative backscatter (RB)
RB = RL + TL;


%% Save into data structure:

% Loop through beams and save
for n = 1:4
    Data.([dataModeWord,'_RelBackscatterBeam',num2str(n)]) = RB(:,:,n);
end

% Also save the average of all beams
Data.([dataModeWord,'_RelBackscatter']) = nanmean(RB,3);

% ...and the attenuation coefficient, alpha
Data.attenuation_coefficient = alpha;

end