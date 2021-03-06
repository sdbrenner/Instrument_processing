function [Data,magOffset] = sigMagCorrection2D( Data, mode, magOffset, makePlots )
% SIGMAGCORRECTION2D "Hard-iron" magenetic compass correction for Nortek
% Signature-series acoustic doppler current profilers
%
%   Data = sigMagCorrection2D( Data ) calculates and applies a "hard-iron"
%   magnetic correction to Signature-series heading data for Average mode.
%   This correction completely replaces the Data.Average_Heading variable
%   instead of just modifying it, so other corrections applied (e.g.
%   declination offset) are wiped out.
%
%   [ Data , magOffset ] = sigMagCorrection2D( Data ) outputs a 1x2 vector
%   magOffset containing the x and y magnetic offset values calculated.
%
%	[ ... ] = sigMagCorrection2D( Data, mode) allows for specification of
%	the input data mode as 'avg', 'burst', or 'ice' (corresponding to
%	Average, Burst, or AverageIce structure variables).  The function can
%	act on multiple data types by including different modes by including a
%	cell array of modes: e.g. {'avg','ice'}
%
%   [ ... ] = sigMagCorrection2D( Data, mode, magOffset) applies the x and
%   y magnetic offset values from the vector 'magOffset' instead of
%   calculting them.
%
%   [ ... ] = sigMagCorrection2D( ..., makePlots ) defines whether or not
%   the function should make plots showing the magnetic circle and
%   corrected v.s. uncorrected heading.  The variable 'makePlots' is
%   specified as a booleen true (1) or false (0). By default the function
%   will make plots (makePlots=1) if it is calculating the magnetic offset
%   values, but will not make plots (makePlots=0) if an offset value is
%   specified.
%
%   Notes:  
%   (1) This function is developed to operate on Data structures that are
%   output by converting raw .ad2cp data to .mat files using MIDAS
%   software.  Data converted with Signature Deployment software may not
%   have matching variable names.
%
%   S.D.Brenner, 2019
%
%   Updated
%   Apr., 2020:     Corrected Earth reference frame horizontal projection
%                   (per algorithm from Dan Torres, WHOI)

%% 

% By default:
%   Calculate a magnetic offset correction
%   Calculate correction from/apply correction to 'Average' mode
%   Make plots showing circle fit
correctionCalculated = 0;
if nargin < 4 || isempty(makePlots); makePlots = 1;       end
if nargin < 2 || isempty(mode);      mode = 'avg';    end

% If a magnetic offset vector is provided, don't perform calculation
% Also, the default plot behaviour in this case is to NOT make plots.
% ( note: to make plots, a magnetic circle radius is required, which can't 
%    be found from the offset values; thus a circle fit must still be
%    performed ).
if nargin >= 3 && ~isempty(magOffset)
    MoX = magOffset(1); 
    MoY = magOffset(2);
    correctionCalculated = 1;
    if nargin < 4 || isempty(makePlots); makePlots = 0; end
else
    magOffset = [];
end


% Parse mode choice
%   ( Note, 'mode' options could have instead been the 'dataWordChoices'
%     values, but instead are 'modeChoices' to be consistent with other
%     Nortek and Signature codes)
modeChoices = {'avg','ice','burst'};
dataWordChoices = {'Average','AverageIce','Burst'};
[lia,locb] = ismember( lower(mode) , modeChoices );
if ~lia
    error('The input variable ''mode'' must be one of: ''avg'', ''ice'', or ''burst''');
elseif length(lia)>1
    % If multiple mode words are entered, recursively run this script for
    % each of the individually (this may break something)
    for n = 1:length(lia)
        modeN = modeChoices{locb(n)};
        [Data,magOffsetN] = sigMagCorrection2D( Data, modeN, magOffset, makePlots );
        magOffset(n,:) = magOffsetN;
    end
    return;
else
    dataModeWord = dataWordChoices{locb};
end


%% Extract data

if isfield( Data, [dataModeWord '_HeadingRaw'] )
    heading = Data.( [dataModeWord '_HeadingRaw'] );
else
    heading = Data.( [dataModeWord '_Heading'] );
end

pitch = Data.( [dataModeWord '_Pitch'] );
roll = Data.( [dataModeWord '_Roll'] );

sigMagX = Data.( [dataModeWord '_MagnetometerX'] );
sigMagY = Data.( [dataModeWord '_MagnetometerY'] );
sigMagZ = Data.( [dataModeWord '_MagnetometerZ'] );

%% Rotate from Signature reference frame to flat (Earth reference frame):
% i.e.: remove pitch/roll

% "Naive" correction (previously applied; does not account for z-component
% of magnetic field):
nMagX = sigMagX./cosd(pitch);
nMagY = sigMagY./cosd(roll);

% "Robust" correction (from Dan Torres)
sinP = sind(pitch);
cosP = cosd(pitch);
sinR = sind(roll);
cosR = cosd(roll);
sinPsinR = sinP.*sinR;
sinPcosR = sinP.*cosR;
magX = sigMagX.*cosP - sigMagY.*sinPsinR - sigMagZ .* sinPcosR;
magY = sigMagY.*cosR - sigMagZ.*sinR;


% Rotate using transformation matrix from M Le Menn et al., 2014 
% ( doi:10.1088/0957-0233/25/8/085801 )
% * in order to match result from Dan Torres, use negative roll
% sP = reshape( sind(pitch), 1,1,[] );
% cP = reshape( cosd(pitch), 1,1,[] );
% sR = reshape( sind(-roll), 1,1,[] );
% cR = reshape( cosd(-roll), 1,1,[] );
% Z = zeros(size(sP));
% TM = [cP, (sR.*sP), (-cR.*sP) ; Z, cR, sR];
% sigMag = [sigMagX.' ; sigMagY.'; sigMagZ.' ];
% for n = 1:size(TM,3)
%     mag(:,n) = TM(:,:,n) * sigMag(:,n) ;
% end
% magX = mag(1,:);
% magY = mag(2,:);


%% Remove outlier data:
% the circle fit seems somewhat sensitive to outliers, so here I am making
% a try at removing them. This is unsophisticated, so it may result in
% other problems.

% magX = filloutliers(magX,'linear','mean');
% magY = filloutliers(magY,'linear','mean');

%% Calculate correction

if ~correctionCalculated
% Fit a circle to the data to find offset (Mo), and corrected magnetic
% circle (Mc) radius 
    [MoX,MoY,McR] = circfit(magX,magY);
    magOffset(1) = MoX;
    magOffset(2) = MoY;
elseif correctionCalculated && makePlots
% To make plots, a magnetic circle radius is required, which can't be found
% from the offset values
    R = sqrt( (magX-MoX).^2 + (magY-MoY).^2 );
    McR = median(R);
end

if makePlots
    % Define coordinates of ideal corrected circle (Mc)
    McTH = linspace(-pi,pi);
    [McX,McY] = pol2cart( McTH, McR);
    MX = McX + MoX;
    MY = McY + MoY;
end



%% Apply correction

% Apply to correction to heading data
correctedHeading = mod( atan2d(magY-MoY,magX-MoX) , 360 );
Data.([dataModeWord  '_HeadingRaw']) = heading;
Data.([dataModeWord  '_Heading']) = correctedHeading;


% Check if the Data structure includes an AHRS rotation matrix:
ahrsEnabled = isfield(Data,[dataModeWord  '_AHRSRotationMatrix']);

% If the machine has an AHRS system, the heading correction needs to be
% incorporated into the AHRS rotation matrix. I will follow a procedure 
% similar to what is available in a Nortek-provided code for declination 
% correction, but replace the declination angle with the correction angle.
% This needs to be done separately for each timestamp.
if ahrsEnabled 
    numTimes = length( Data.([dataModeWord,'_Heading']));
    for n = 1:numTimes
        % extract and organize AHRS rotation matrix
        rotMatAHRS = ...
 transpose(reshape(Data.([dataModeWord  '_AHRSRotationMatrix'])(n,:),3,3));
        % calculate correction angle
        corAng = mod( correctedHeading(n) - heading(n), 360);
        % construct and appy correction rotation matrix
        corMat = [  cosd(corAng) , sind(corAng) , 0 ;
                   -sind(corAng) , cosd(corAng) , 0 ;
                         0,             0,        1 ];
        corRotMatAHRS = corMat *  rotMatAHRS;            
        % update AHRS rotation matrix
        Data.([dataModeWord  '_AHRSRotationMatrix'])(n,:) = ... 
                                     reshape(transpose(corRotMatAHRS),1,9);
    end
end





%% Plot

if makePlots
    figure; clf;
    subplot(1,2,1)
    hold on;
    % plot magnetometer data
    scatter( magX,magY,4, 'k','filled','MarkerFaceAlpha',0.25);
    % Add circle fit
    plot( McX + MoX , McY + MoY , 'r','linewidth',2);
    % Add center marking and offset
    plot( [0,MoX],[0,MoY],'bo-');
    xline( MoX,'r');
    yline( MoY,'r');
    % Set axis and labels
    set(gca,'xlim',[MoX-1.5*McR,MoX+1.5*McR],...
            'ylim',[MoY-1.5*McR,MoY+1.5*McR] );
    daspect([1,1,1]);
    grid on;
    xlabel('Mag. X (rotated)');
    ylabel('Mag. Y (rotated)');
    text(0.05,0.95,'(a)','units','normalized','fontsize',12);
    
    subplot(1,2,2);
    hold on;
%     one2one = [0,360];
%     scatter(heading, correctedHeading , 4 ,...
%             'filled','MarkerFaceAlpha',0.25);
%     plot( one2one, one2one, 'k--','linewidth',2);
    scatter(heading, mod(correctedHeading-heading+180,360)-180 , 4 ,'k',...
            'filled','MarkerFaceAlpha',0.25);    
    yline( 0 , 'k--','linewidth',2);
    set(gca,'xlim',[0,360],'xtick',0:90:360);
    set(gca,'ylim',[-30,30]);
%     set(gca,'ylim',[0,360],'ytick',0:90:360 );
%     daspect([1,1,1]);
    grid on;
    xlabel('Reported heading [deg]');
    ylabel('Heading correction [deg]');
    text(0.05,0.95,'(b)','units','normalized','fontsize',12);
end

end
%% Helper functions

function   [xc,yc,R,a] = circfit(x,y)
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%  By:  Izhak Bucher 25/oct /1991, 

    x=x(:); y=y(:);
   a=[x y ones(size(x))]\[-(x.^2+y.^2)];
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
   
end



