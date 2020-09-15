function fig = showPCAreconstruction(D,k,npc)
%SHOWPCARECONSTRUCTION Make bar graph of channelwise error from PCA reconstruction
%
%  fig = kal.showPCAreconstruction(D,k,npc);
%
% Inputs
%  D     - Table
%  k     - Index into table row
%  npc   - (default is all PCs) | scalar (number of PCs to use in
%           reconstruction of original data)
%
% Output
%  fig   - Figure handle to bar graph of error estimates
%
% See also: kal, analyze, analyze.jPCA, analyze.dynamics

% parse input
if nargin < 2
   k = 116; % Index of recording (116 = RC-43:Grasp:D26, high R2_best)
end

% Get scores and coefficients for this recording
X = D.Summary{k}.PCA.scores; % RC-43:Grasp:D26 PCA scores (states)
H = D.Summary{k}.PCA.vectors_all; % RC-43:Grasp:D26 PCA coefficients
N_PC = size(H,2);
if nargin < 3
   npc = N_PC;
end

% Get means as well
mu = D.Summary{k}.PCA.mu; % RC-43:Grasp:D26 PCA means

% Get dataset spike rates
Z = vertcat(D.Projection{k}.data); % RC-43:Grasp:D26 spike rates

% Cross-trial mean
n = numel(D.Projection{k});
mu_z = repmat(nanmean(cat(3,D.Projection{k}.data),3),n,1);
% sd_z = repmat(nanstd(cat(3,D.Projection{k}.data),[],3),n,1);

% Do subtraction
% Z0 = ((Z - mu_z)./sd_z)'; % "Observations"
Z0 = (Z - mu_z).';

% data = scores * coeff' + mu; 
% == 
% data' = coeff * scores' + mu';
% ===
% Z = H * X + v;

Zhat = H(:,1:npc) * X(:,1:npc)' + mu'; 
err = nanmean((Z0 - Zhat).^2,2);
err_max = nanmax(err);
if err_max < 1e-20
   YL = [eps^2 1];
   YS = 'log';
elseif err_max < 1e-10
   YL = [1e-15 1]; 
   YS = 'log';
elseif err_max < 1e-5
   YL = [1e-10 1];
   YS = 'log';
elseif err_max < 1e-3
   YL = [1e-5 1];
   YS = 'log';
elseif err_max < 1e-1
   YL = [1e-3 1];
   YS = 'log';
elseif err_max < 1
   YL = [0 1];
   YS = 'linear';
else
   YL = [0 1].*2.*nanmean(nanvar(X));
   YS = 'linear';
end
   
   
fig = figure('Name','Proof of Concept: PCA fit','Color','w');
ax = axes(fig,'NextPlot','add',...
   'XLim',[0 (numel(err)+1)],...
   'YLim',YL,...
   'YScale',YS,...
   'XColor','k','YColor','k',...
   'LineWidth',1.5,'FontName','Arial');
bar(ax,1:numel(err),err,1,...
   'EdgeColor','none',...
   'FaceColor','k',...
   'DisplayName','Reconstruction Error');
legend(ax,...
   'TextColor','black',...
   'Location','Northeast',...
   'FontName','TimesNewRoman');
ylabel(ax,'Mean Square Error','FontName','Arial','Color','k');
xlabel(ax,'Channel','FontName','Arial','Color','k');
title(ax,sprintf('PCA Reconstruction (%d / %d PCs)',npc,N_PC),...
   'FontName','Arial','Color','k');


end