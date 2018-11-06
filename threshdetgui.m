function threshdetgui
% ** function threshdetgui
% Graphical user interface for the detection of neuronal 'events' (e.g.
% action potentials) and 'bursts' (e.g. Up States) in time series data via
% threshold methods.

% -------------------------------------------------------------------------
% Version 5.10.1, November 2018
% (C) Harald Hentschke (University Hospital of Tuebingen)
% -------------------------------------------------------------------------

global sp wp

% Here, layout of all gui elements (knobs, buttons, etc.) should be done
labelscale('fontSz',8,'scaleFac',1.0,'lineW',.25,'markSz',6); 
% standard dimension of edit fields
editw=.08;
edith=.032;
% standard dimensions of buttons
butt1x=.08;
butt1y=.04;
butt2x=.16;
butt2y=.08;

% elementary delta y given height of text fields and smaller buttons
dy=.036;
% standard text field width (at least as wide as big buttons)
textw=.16;
% standard margin for major divisions of figure (like plots and groups of buttons)
smarg=.025;

% left alignment lines in ascending order
la1=smarg;
% left border left-aligned subplots
la2=la1+textw+smarg;
% left border of right subplot 
la3=(1-la2+smarg)/2+la2;

% vertical separators from top to bottom
ba1=.96;
ba2=.55;
ba11=.33;

% font sizes
fsz=8;
fsz_big=11;

% if gui exists from previous session, delete it
tmph=findobj('Tag','threshdetgui','type','figure');
if ~isempty(tmph)
  delete(tmph);
end

% the collection of routines associated with most button callbacks
funcH1=@threshdetguifunc;

% create GUI window
H0 = figure('Units','normalized', ...
   'Color',[.9 .9 .9], ...
   'Name','Event Detection via Threshold', ...
   'NumberTitle','off', ...
   'Position',[0.005 0.25 0.99 0.7], ...
   'Tag','threshdetgui',...
   'MenuBar','none',...
   'Toolbar','figure',...
   'DeleteFcn',{funcH1,{'done'}}...
);

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-1*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','read data channel', ...
  'TooltipString','',...
  'Tag','readDataBttn',...
  'callback', {funcH1,{'readData'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-2*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','options', ...
  'TooltipString','',...
  'Tag','optionsBttn', ...
  'callback', {funcH1,{'openOptionsDialog'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-3*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','data preconditioning', ...
  'TooltipString','',...
  'Tag','precondBttn', ...
  'callback', {funcH1,{'preconditionRawData'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-4*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','detect events', ...
  'ForegroundColor',[0.4 0.4 1.0], ...
  'TooltipString','',...
  'Tag','detEvtBttn', ...
  'callback', {funcH1,{'detEvt'}});  


uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-5*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','event amp', ...
  'TooltipString','computation of PSC amplitudes',...
  'Tag','evtAmpBttn', ...
  'callback', {funcH1,{'detEvtAmp'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-6*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','separate pops (PCA)', ...
  'TooltipString','',...
  'Tag','sepEvtBttn', ...
  'callback', {funcH1,{'sepEvt'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-7*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','inspect/purge pops', ...
  'TooltipString','allows you to either delete separate(d) populations of events or delete single events',...
  'Tag','treatEvtBttn', ...
  'callback', {funcH1,{'treatEvt'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-8*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','detect bursts', ...
  'ForegroundColor',[.9 0.3 0.3], ...
  'TooltipString','',...
  'Tag','detBurstBttn', ...
  'callback', {funcH1,{'detBurst'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-9*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'foregroundcolor','m',...
  'style','pushbutton',...
  'String','experimental', ...
  'TooltipString','',...
  'Tag','saveResBttn', ...
  'callback', {funcH1,{'experimental'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-10*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','save results', ...
  'TooltipString','',...
  'Tag','saveResBttn', ...
  'callback', {funcH1,{'saveResults'}});  

% this is a hidden, inactive button that will become visible and enabled
% only in batch mode, but whose status will be enquired with each loop run
% of threshdetguifunc
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-11*butt2y butt2x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'foregroundcolor','y',...
  'backgroundcolor','k',...
  'style','togglebutton',...
  'String','cancel batch job', ...
  'TooltipString','',...
  'visible','off',...
  'enable','off',...
  'Tag','cancelBatchJobBttn');

  
% Set up the subplots and don't forget to tag them so they can be identified 
% (plots need larger margins because of axes labels)
% - overview raw data
sp.rawOv.axH=subplot('position',[la2  0.75+smarg  1-la2-smarg  0.25-2*smarg]);
set(sp.rawOv.axH,'tag','rawOverview','NextPlot','add');
% - excerpt raw data
sp.rawExc.axH=subplot('position',[la2  0.5+smarg  1-la2-smarg  0.25-2*smarg]);
set(sp.rawExc.axH,'tag','rawExcerpt');
% - cutouz
sp.cutout.axH=subplot('position',[la2  3*smarg  0.5-la2-smarg  0.5-5*smarg]);
set(sp.cutout.axH,'tag','cutouts');
% - iei
sp.iei.axH=subplot('position',[la3  3*smarg  0.5-la2-smarg  0.5-5*smarg]);
set(sp.iei.axH,'tag','iei');

% field .isBatchMode has to be set before calling threshdetguifunc because
% in deployed mode one of the first code lines in threshdetguifunc requires
% its existence despite the short-circuit &&
wp.isBatchMode=false;

% threshdetguifunc(H0,[],{'init','openOptionsDialog'},'sp',sp);
threshdetguifunc(H0,[],{'init','openOptionsDialog'});

