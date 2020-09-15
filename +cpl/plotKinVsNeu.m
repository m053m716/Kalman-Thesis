function fig = plotKinVsNeu(snips,tNeu,worldPoints,tKin,varargin)
%% PLOTKINVSNEU   Plot kinematic movement vs world point data
%  
%  fig = PLOTKINVSNEU(snips,tSnips,worldPoints,tWorld);
%  fig = PLOTKINVSNEU(snips,tSnips,worldPoints,tWorld,'NAME',value,...);
%
%  --------
%   INPUTS
%  --------
%    snips     :     Table with snippets of neural data aligned to grasps.
%
%   tNeu       :     Times corresponding to neu snippets (relative ot grasp
%                       alignment)
%
%  worldPoints :     Kinematic data for all markerless data, aligned to
%                       grasps.
%
%   tKin       :     Timing for each point in worldPoints series (at a
%                       different rate that neural).
%
%  varargin    :     (Optional) 'NAME', value input argument pairs
%
%  --------
%   OUTPUT
%  --------
%     fig      :     Handle to figure that contains subplots comparing the
%                       kinematic trajectories to neural trajectories for
%                       grasps.
%
% By: Max Murphy  v1.0  11/21/2018  Original version (R2017b)

%% DEFAULTS
VAR_NAMES = {'d1_d';'d1_p';'d2_d';'d2_p'};
VAR_IDX = [1,3,4,6];

KIN_DIM = {'IC_1';'IC_2';'IC_3'};

NEU = 'SG';
NEU_YLIM = [-5.0 5.0];
NEU_CHAN = 16;

N_LINE = 3;
REACH_IDX = nan;
DO_ICA = true;
LINE_COL = repmat(linspace(0,0.75,N_LINE).',1,3);
LINE_W = fliplr(linspace(1,2.5,N_LINE));

AUTO_LIMS = false;

MOVIE = nan;

%% PARSE VARARGIN
for iV = 1:2:numel(varargin)
   eval([upper(varargin{iV}) '=varargin{iV+1};']);
end

%% 
fig = figure('Name','Kinematics vs Neural',...
       'Units','Normalized',...
       'Color','w',...
       'Position',[0.2 0.2 0.6 0.7]);
    
% Determine number of subplots
nRow = 7;
nCol = numel(VAR_IDX);

% Define neural axis first
neuAx = subplot(nRow,nCol,1:(2*nCol));
neuAx.NextPlot = 'add';
neuAx.XLimMode = 'manual';
neuAx.XTick = [];
neuAx.XLim = [min(tNeu) max(tNeu)];
neuAx.YLimMode = 'manual';
neuAx.YLim = NEU_YLIM;

% xlabel('Time (ms)','FontName','Arial','Color','k');
title(sprintf('%s Neural Signal (%s-%s)',...
   NEU,snips.Probe{NEU_CHAN},snips.Channel{NEU_CHAN}),...
   'FontName','Arial','FontSize',14,'Color','k');
ylabel('Z-score','FontName','Arial','Color','k');

% Then define kinematic axes
if AUTO_LIMS
   [XL,YL,ZL] = cpl.getWPaxesLims(worldPoints);
%    XL = [5 35];
%    ZL = [65 110];
else
   XL = [-4 4];
   YL = [-5 5];
   ZL = [-4 4];
end

LIMS = [XL;YL;ZL];
kinAx = cell(4,numel(VAR_IDX));

for iC = 1:3
   kinAx{iC,1} = subplot(nRow,nCol,((iC+1)*nCol+1):((iC+2)*nCol));
   kinAx{iC,1}.NextPlot = 'add';
   kinAx{iC,1}.XLimMode = 'manual';
   kinAx{iC,1}.XLim = [min(tKin) max(tKin)];
   kinAx{iC,1}.YLimMode = 'manual';
   kinAx{iC,1}.YLim = LIMS(iC,:);
   
   if iC==1
      title('Independent Components','FontName','Arial','Color','k',...
         'FontSize',16);
   end

   if iC==3
      xlabel('Time (ms)','FontName','Arial','Color','k');
   end


   ylabel(KIN_DIM{iC},'FontName','Arial','Color','k');


end


% Last, define special 2D axes
for iA = 1:2:nCol
   kinAx{4,iA} = subplot(nRow,nCol,[(nCol*5+iA) (nCol*6+iA)]);
   kinAx{4,iA}.NextPlot = 'add';
   
   kinAx{4,iA}.XLimMode = 'manual';
   kinAx{4,iA}.YLimMode = 'manual';
%    kinAx{4,iA}.XLim = XL;
   kinAx{4,iA}.XLim = [5 35];
%    kinAx{4,iA}.YLim = ZL;
   kinAx{4,iA}.YLim = [65 110];
   kinAx{4,iA}.YDir = 'reverse';
   
   xlabel('X-Dim (mm)','FontName','Arial','Color','k');
   if iA==1
      ylabel('Z-Dim (mm)','FontName','Arial','Color','k');
   end
   
end
    
for iA = 2
   kinAx{4,iA} = subplot(nRow,nCol,[(nCol*5+iA) (nCol*6+iA)]);
   kinAx{4,iA}.NextPlot = 'add';
   
   kinAx{4,iA}.XLimMode = 'manual';
   kinAx{4,iA}.YLimMode = 'manual';
   kinAx{4,iA}.XLim = XL;
%    kinAx{4,iA}.XLim = [5 35];
   kinAx{4,iA}.YLim = YL;
%    kinAx{4,iA}.YLim = [65 110];
%    kinAx{4,iA}.YDir = 'reverse';
   
   xlabel('IC_1','FontName','Arial','Color','k');
   ylabel('IC_2','FontName','Arial','Color','k');

   
end

for iA = 4
   kinAx{4,iA} = subplot(nRow,nCol,[(nCol*5+iA) (nCol*6+iA)]);
   kinAx{4,iA}.NextPlot = 'add';
   
   kinAx{4,iA}.XLimMode = 'manual';
   kinAx{4,iA}.YLimMode = 'manual';
   kinAx{4,iA}.XLim = XL;
   kinAx{4,iA}.YLim = ZL;
   
   xlabel('IC_1','FontName','Arial','Color','k');
   ylabel('IC_3','FontName','Arial','Color','k');

end

nTrial = size(worldPoints,1);

neuLine = [];
kinLine = cell(3,1);
kinMark = cell(1,numel(VAR_IDX));
worldZ = cpl.getKinICs(worldPoints,'VAR_IDX',VAR_IDX,'DO_ICA',DO_ICA);


if ~isnan(MOVIE)
   V = VideoWriter(MOVIE,'Motion JPEG AVI');
   V.FrameRate = 0.5;
   open(V);
end

if isnan(REACH_IDX)
   trialVec = 1:nTrial;
else
   trialVec = reshape(REACH_IDX,1,numel(REACH_IDX));
end

for iTrial = trialVec
   plot(neuAx,tNeu,snips.(NEU){NEU_CHAN}(iTrial,:).',...
      'Color',LINE_COL(1,:),...
      'LineWidth',LINE_W(1),'UserData',true);
   c = get(neuAx,'Children');
   updateLineCols(c,LINE_COL,LINE_W,N_LINE);
   
   [~,idx] = min(snips.(NEU){NEU_CHAN}(iTrial,:));
   tMin = tNeu(idx);
   
   
   if isempty(neuLine)
      neuLine = line(neuAx,[tMin tMin],NEU_YLIM,'Color','m','LineWidth',2,'Marker','v',...
         'MarkerFaceColor','m','MarkerIndices',1,'UserData',false);
   else
      neuLine.XData = [tMin tMin];
   end
   

   for iC = 1:3
      plot(kinAx{iC,1},tKin,worldZ{iTrial}(:,iC),...
         'LineWidth',LINE_W(1),'Color',LINE_COL(1,:),'UserData',true);
      c = get(kinAx{iC,1},'Children');
      updateLineCols(c,LINE_COL,LINE_W,N_LINE);
      if isempty(kinLine{iC,1})
         kinLine{iC,1} = line(kinAx{iC,1},[tMin tMin],LIMS(iC,:),'Color','m',...
            'LineWidth',2,'LineStyle',':','UserData',false);
      else
         kinLine{iC,1}.XData = [tMin tMin];
      end

   end
      
   for iA = 1:2:nCol
      plot(kinAx{4,iA},...
           worldPoints{iTrial}{VAR_IDX(iA)}(:,1), ...
           worldPoints{iTrial}{VAR_IDX(iA)}(:,3),'LineWidth',LINE_W(1),...
           'Color',LINE_COL(1,:),'UserData',true);
      c = get(kinAx{4,iA},'Children');
      [~,idx] = min(abs(tKin - tMin));
      if isempty(kinMark{iA})
         kinMark{iA} = plot(kinAx{4,iA},worldPoints{iTrial}{VAR_IDX(iA)}(idx(1),1),...
            worldPoints{iTrial}{VAR_IDX(iA)}(idx(1),3),...
            'LineStyle','none','Marker','o',...
            'MarkerFaceColor','m','MarkerEdgeColor','r',...
            'MarkerSize',10,'UserData',false);
      else
         kinMark{iA}.XData = worldPoints{iTrial}{VAR_IDX(iA)}(idx(1),1);
         kinMark{iA}.YData = worldPoints{iTrial}{VAR_IDX(iA)}(idx(1),3);
      end
      
      updateLineCols(c,LINE_COL,LINE_W,N_LINE);

      
   end
   
   for iA = 2
      plot(kinAx{4,iA},...
           worldZ{iTrial}(:,1), ...
           worldZ{iTrial}(:,2),'LineWidth',LINE_W(1),...
           'Color',LINE_COL(1,:),'UserData',true);
      c = get(kinAx{4,iA},'Children');
      [~,idx] = min(abs(tKin - tMin));
      if isempty(kinMark{iA})
         kinMark{iA} = plot(kinAx{4,iA},worldZ{iTrial}(idx(1),1),...
            worldZ{iTrial}(idx(1),2),...
            'LineStyle','none','Marker','o',...
            'MarkerFaceColor','m','MarkerEdgeColor','r',...
            'MarkerSize',10,'UserData',false);
      else
         kinMark{iA}.XData = worldZ{iTrial}(idx(1),1);
         kinMark{iA}.YData = worldZ{iTrial}(idx(1),2);
      end
      
      updateLineCols(c,LINE_COL,LINE_W,N_LINE);  
   end
   
   for iA = 4
      plot(kinAx{4,iA},...
           worldZ{iTrial}(:,1), ...
           worldZ{iTrial}(:,3),'LineWidth',LINE_W(1),...
           'Color',LINE_COL(1,:),'UserData',true);
      c = get(kinAx{4,iA},'Children');
      [~,idx] = min(abs(tKin - tMin));
      if isempty(kinMark{iA})
         kinMark{iA} = plot(kinAx{4,iA},worldZ{iTrial}(idx(1),1),...
            worldZ{iTrial}(idx(1),3),...
            'LineStyle','none','Marker','o',...
            'MarkerFaceColor','m','MarkerEdgeColor','r',...
            'MarkerSize',10,'UserData',false);
      else
         kinMark{iA}.XData = worldZ{iTrial}(idx(1),1);
         kinMark{iA}.YData = worldZ{iTrial}(idx(1),3);
      end
      
      updateLineCols(c,LINE_COL,LINE_W,N_LINE);  
   end
   
   if ~isnan(MOVIE)
      I = getframe(fig);
      writeVideo(V,I);
   end
   
   
   
end

if ~isnan(MOVIE)
   close(V);
end


   function updateLineCols(c,lc,lw,nl)
      keep_idx = true(numel(c),1);
      for ii = 1:numel(c)
         keep_idx(ii) = c(ii).UserData;
      end
      cc = c(keep_idx); % remove annotation lines      
      
      if (numel(cc)>nl)
         delete(cc(end));
         cc(end) = [];
      end
      
      for iLine = 1:numel(cc)
         cc(iLine).LineWidth = lw(iLine);
         cc(iLine).Color = lc(iLine,:);
      end
      
   end


end