function H0=cutdealgui
% ** function H0=cutdealgui
% Graphical user interface for dealing with 'cutouts' of events. It is part
% of threshdetgui and not self-sufficient. fh is the handle to the GUI
% figure.

% Here, layout of all gui elements (buttons, etc.) should be done

% standard dimension of edit fields
editw=.08;
edith=.032;
% standard dimensions of buttons
butt1x=.12;
butt1y=.05;
butt2x=.16;
butt2y=.10;

% elementary delta y given height of text fields and smaller buttons
dy=.036;
% standard text field width (at least as wide as big buttons)
textw=.16;
% standard margin for major divisions of figure (like plots and groups of buttons)
smarg=.025;

% left alignment lines in ascending order
la1=smarg;
% left border left-aligned subplots
la2=la1+butt1x+smarg;
% left border of right subplot 
la3=(1-la2+smarg)/2+la2;

% vertical separators from top to bottom
ba1=.96;
ba2=.55;
ba11=.33;

% font sizes
fsz=8;
fsz_big=10;

% if gui exists from previous session, delete it
tmph=findobj('name', 'Processing of Cutouts', 'type', 'figure');
if ~isempty(tmph)
  delete(tmph);
end

% the collection of routines associated with most button callbacks
funcH1=@cutdealguifunc;

% create GUI window
H0 = figure('Units','normalized', ...
   'Color',[0.8 0.8 0.8], ...
   'Name','Processing of Cutouts', ...
   'NumberTitle','off', ...
   'Position',[0.005 0.05 0.99 0.9], ...
   'Tag','CutDeal',...
   'MenuBar','none',...
   'Toolbar','figure'...
);

% Set up the subplots and don't forget to tag them so they can be identified 
% (plots need larger margins because of axes labels)
% - overview 
sp.parall.axH=subplot('position',[la1  0.12+smarg  0.6  1-0.12-2*smarg]);
set(sp.parall.axH,'tag','parallelPlot','NextPlot','add');
% - overlay, mean + fits
sp.overlay.axH=subplot('position',[la1+0.6+smarg  0.6  1-la1-0.6-2*smarg  0.4-smarg]);
set(sp.overlay.axH,'tag','overlay');


% ---- parallel plot control buttons
h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-18*butt1y 0.3 butt1y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','<< previous <<', ...
  'TooltipString','',...
  'Tag','prevBttn', ...
  'callback', {funcH1,{'plotPrev'}});  

h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+0.3 ba1-18*butt1y 0.3 butt1y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','>> next >>', ...
  'TooltipString','',...
  'Tag','nextBttn', ...
  'callback', {funcH1,{'plotNext'}});  

h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la1+0*butt1x ba1-19*butt1y butt1x*.5 butt1y], ... 
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'String','rows cols', ...
  'Tag','setPlotLayout',...
  'callback', {funcH1,{'setPlotLayout'}});  

h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+1*butt1x ba1-19*butt1y butt1x butt1y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','deflate', ...
  'TooltipString','increase vertical offset between cutouts',...
  'Tag','incYOffsBttn', ...
  'callback', {funcH1,{'incYOffs'}});  

h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+2*butt1x ba1-19*butt1y butt1x butt1y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','pump up', ...
  'TooltipString','decrease vertical offset between cutouts',...
  'Tag','decYOffsBttn', ...
  'callback', {funcH1,{'decYOffs'}});  



% ---- analysis/action buttons
h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+0.6+smarg 0.6-2*butt1y butt1x butt1y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','average', ...
  'TooltipString','average of cutouts will be displayed upon button press',...
  'Tag','averageBttn', ...
  'callback', {funcH1,{'averageEvt'}});  

h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+0.6+smarg 0.6-9*butt1y butt1x butt1y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','remove tagged', ...
  'TooltipString','cutouts tagged as bad will be removed upon button press',...
  'Tag','rmBttn', ...
  'callback', {funcH1,{'rmEvt'}});  

% ----- OK button
h1 = uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[.85 ba1-18*butt1y butt1x butt2y], ...
  'FontSize',fsz_big, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','DONE', ...
  'TooltipString','And now for something completely different...',...
  'ListboxTop',0, ...
  'Tag','doneBttn', ...
  'callback', {funcH1,{'done'}});  


cutdealguifunc(H0,[],{'init'},'sph',sp);
% cutdealguifunc(H0,[],{'init','openOptionsDialog'},'sph',sp);
