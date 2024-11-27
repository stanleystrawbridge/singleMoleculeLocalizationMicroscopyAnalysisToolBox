
% INSTRUCTIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Endesfelder, U., Malkusch, S., Fricke, F. et al. A simple method to
% estimate the average localization precision of a single-molecule
% localization microscopy experiment. Histochem Cell Biol 141, 629–638
% (2014). https://doi.org/10.1007/s00418-014-1192-3
%
% Assumption: Localizations are normally distributed (localization
% positions are gaussian)
%
% Limitations: Cannot be used for, 
% - dense targets where single-molecule distributions are not resolved
%
%- Experiments using fluorophores that are not localized many times (e.g.,
%   irreversibly switching fluorophores used for photoactivated localization
%   microscopy (PALM))
%
% OR 
%
% - or time-critical experiments where, e.g., fast dynamics of living cells
%   are involved, hinder accurate determination of σSMLM
%
% - This Gaussian approximation is valid for isotropic symmetric
%   single-molecule signals, when effects like light scattering in thicker
%   samples, a large refraction index mismatch, dipole orientations or
%   diffraction rings due to, e.g., larger imaging depths are excluded
%
% Nearest Neighbors in adjacent frames NN_{adjf}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==========================================================================
% USER INPUT BEGIN
%==========================================================================
% Uses the '*_trackPositions.csv' from trajectory analysis script

inputPath = ['G:\Shared drives\Srinjan LAB\Theses_Publications_in_progress\Publication_Shuchi_HOXA9\smlm_analysis\20230630_localizations\raw\100ms_out']


% Is the data 'live' or 'fixed'?
imageType = 'fixed'; %'live' WORK IN PROGRESS, DO NOT CHANGE

%==========================================================================
% USER INPUT END
%==========================================================================

obj = calculateSmlmPrecision(inputPath,imageType);







