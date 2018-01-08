function threshdet_optgui
% ** function threshdet_optgui
% gui for setting analysis & general options within threshdetgui

% ***************** WATCH OUT *********************************************
% whenever you (as the programmer) make any change to this gui, including 
% and notably the callbacks of the buttons, saved parameter files must be 
% recreated by the user because these parameter files are matlab
% 'fig'-files which store everything, including the callbacks. 'Recreating'
% does not necessarily mean retyping all values from scratch. In case only
% the callbacks changed, do the following:
% open options dialogue
% - load old parameter file
% - hit 'apply' button
% - close window manually if it does not close automatically
% - reopen options dialogue
% - save (same or other name), done!
% *************************************************************************

% standard dimension of edit fields
editw=.12;
edith=.036;

% standard dimensions of buttons
butt1x=.08;
butt1y=.036;
butt2x=.16;
butt2y=.072;

% elementary delta y given height of text fields and smaller buttons
dy=.038;
% standard text field width 
textw=.2;
% left & right alignment lines in ascending order
la1=.01;
la2=la1+textw*1.05;
la3=la2+editw+.01;
la4=la3+textw*1.05;
la5=la4+editw+.01;
% vertical separators from top to bottom
ba1=.96;
ba2=.55;
ba11=.33;

% font sizes
fsz=9;
fsz_big=13;

% callback for buttons
funcH=@threshdetguifunc;

% if gui exists from previous session, delete it
tmph=findobj('name', 'options', 'type', 'figure');
if ~isempty(tmph)
  delete(tmph);
end


% create (new) GUI window
H0 = figure('Units','normalized', ...
   'Color',[0.8 0.8 0.8], ...
   'Name','options', ...
   'NumberTitle','off', ...
   'Position',[0.05 0.1 0.7 0.7], ...
   'Tag','options',...
   'MenuBar','none',...
   'Toolbar','none'...
);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~ signal preconditioning section
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-0*dy la2-la1+editw edith], ...
  'FontSize',fsz, ...
  'FontWeight','bold', ...
  'BackgroundColor',[0.6 0.6 0.6], ...
  'Style','text', ...
  'String','~~~~~ signal preconditioning', ...
  'Tag','headline1');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-1*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','artifact elimination', ...
  'Tag','text_afElimination');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-1*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'TooltipString','list optional input arguments into function elim_artefact in curly braces (must evaluate to a cell array)',...  
  'Tag','afElimination');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-2*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','lowpass cutoff freq (Hz)', ...
  'Tag','text_loCFreq');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-2*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','loCFreq');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-3*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','highpass cutoff freq (Hz)', ...
  'Tag','text_hiCFreq');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-3*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','hiCFreq');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-4*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','notch filter freq (Hz)', ...
  'Tag','text_notchCFreq');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-4*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'TooltipString','enter the frequency of line hum; the signal will be filtered with a steep bandstop filter; multiple entries like 50 100 150 to get rid of harmonics are possible',...
  'Tag','notchCFreq');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-5*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','differentiator filter', ...
  'Tag','text_differfi');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-5*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'TooltipString','list input arguments, except the first two, into function differfi in curly braces (must evaluate to a cell array)',...  
  'Tag','differfi');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-6*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','downsampling factor', ...
  'Tag','text_sampFac');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-6*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'TooltipString','this is the factor by which the raw data will be downsampled after filtering; e.g. the value of 5 would mean that the first, 6th, 11th, etc. point would be kept',...
  'Tag','sampFac');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-7*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','custom command', ...
  'Tag','text_customCommand');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-8*dy editw 2*edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'max',3,...
  'min',1,...
  'TooltipString','Enter matlab commands you wish to run on raw data after preconditioning above. The variable name is ''d'' (example: d=abs(d);)',...
  'Tag','customCommand');

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-9*dy textw edith], ...
  'FontSize',fsz, ...
  'style','checkbox',...
  'String','cutouts from partly conditioned trace', ...
  'value',0,...
  'TooltipString','check to get cutouts from partly conditioned trace (=lopass filtered & downsampled, but nothing else. By default, conditioned trace will be used)',...
  'Tag','doUncondCutout');  

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~ threshold section
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-10*dy la2-la1+editw edith], ...
  'FontSize',fsz, ...
  'FontWeight','bold', ...
  'BackgroundColor',[0.8 0.8 .2], ...
  'Style','text', ...
  'String','~~~~~ threshold', ...
  'Tag','headline2');
uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-11*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','threshold (absolute)', ...
  'Tag','text_thresh');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-11*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','thresh');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-12*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','threshold (relative)', ...
  'Tag','text_relativeThresh');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-12*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','relativeThresh');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-13*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','burst det: thresh fraction', ...
  'Tag','text_thresh');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-13*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','threshFract');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~ event detection section
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-14*dy la2-la1+editw edith], ...
  'FontSize',fsz, ...
  'FontWeight','bold', ...
  'BackgroundColor',[0.4 0.4 1.0], ...
  'Style','text', ...
  'String','~~~~~ event detection', ...
  'Tag','headline2');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-15*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','detection mode (cross or peak)', ...
  'Tag','text_evtDetMode');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-15*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','evtDetMode');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-16*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','dead time (ms)', ...
  'Tag','text_evDeadT');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-16*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','evDeadT');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-17*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','cutout win (ms)', ...
  'Tag','text_winEvtCutout');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-17*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','winEvtCutout');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~ burst detection section
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-18*dy la2-la1+editw edith], ...
  'FontSize',fsz, ...
  'FontWeight','bold', ...
  'BackgroundColor',[.9 0.3 0.3], ...
  'Style','text', ...
  'String','~~~~~ burst detection', ...
  'Tag','headline2');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-19*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','minimal event width (ms)', ...
  'TooltipString','this parameter specifies the minimal width of an EVENT to be accepted as part of a burst',...  
  'Tag','text_minEventWidth');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-19*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','minEventWidth');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-20*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','maximal gap width (ms)', ...
  'Tag','text_maxBurstGapWidth');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-20*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','maxBurstGapWidth');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-21*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','thresh lag (ms)', ...
  'TooltipString','small signals detected by minor threshold will be ignored if occurring within threshlag ms before signals detected by major threshold',...  
  'Tag','text_threshLag');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-21*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','threshLag');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la1 ba1-22*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','cutout win extension (ms)', ...
  'Tag','text_winBuCutout');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la2 ba1-22*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','winBuCutout');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~ display options
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-0*dy la4-la3+editw edith], ...
  'FontSize',fsz, ...
  'FontWeight','bold', ...
  'BackgroundColor',[0.6 0.6 0.6], ...
  'Style','text', ...
  'String','~~~~~ display options', ...
  'Tag','headline3');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-1*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','length excerpt (ms)', ...
  'Tag','text_excLen');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la4 ba1-1*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','excLen');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-2*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','y limits excerpt (mV)', ...
  'Tag','text_excYLim');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la4 ba1-2*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','excYLim');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-3*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','max # cutouts to plot', ...
  'Tag','text_maxNPlotCutout');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la4 ba1-3*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','maxNPlotCutout');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-4*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','IEI: bin width (ms)', ...
  'Tag','text_ieiBinWid');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la4 ba1-4*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','ieiBinWid');

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-5*dy textw edith], ...  
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','IEI: upper bin limit (ms)', ...
  'Tag','text_ieiBinULimd');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la4 ba1-5*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','ieiBinULim');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~ saving options
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-8*dy la4-la3+editw edith], ...
  'FontSize',fsz, ...
  'FontWeight','bold', ...
  'BackgroundColor',[0.6 0.6 0.6], ...
  'Style','text', ...
  'String','~~~~~ saving options', ...
  'Tag','headline4');

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la3 ba1-9*dy textw edith], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','checkbox',...
  'String','EVENT time stamps', ...
  'value',0,...
  'ForegroundColor',[0.4 0.4 1.0], ...
  'TooltipString','check to have time stamps of events (spikes) saved',...
  'Tag','saveEvt');  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la3 ba1-10*dy textw edith], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','checkbox',...
  'String','EVENT cutouts', ...
  'value',0,...
  'ForegroundColor',[0.4 0.4 1.0], ...
  'TooltipString','check to have event cutouts (fixed length) saved',...
  'Tag','saveEvtCutout');  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la3 ba1-11*dy textw edith], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','checkbox',...
  'String','BURST time stamps', ...
  'value',0,...
  'ForegroundColor',[.9 0.3 0.3], ...
  'TooltipString','check to have timing information (start time, duration) of bursts saved',...
  'Tag','saveBu');  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la3 ba1-12*dy textw edith], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','checkbox',...
  'String','BURST cutouts', ...
  'value',0,...
  'ForegroundColor',[.9 0.3 0.3], ...
  'TooltipString','check to have burst cutouts (variable length) saved',...
  'Tag','saveBuCutout');  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la3 ba1-13*dy textw edith], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','checkbox',...
  'String','preconditioned trace', ...
  'value',0,...
  'TooltipString','check to have preconditioned data trace saved',...
  'Tag','savePrecondData');  

uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-14*dy textw edith], ...
  'FontSize',fsz, ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...  
  'String','string in results file names', ...
  'Tag','text_resFnString');
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la4 ba1-14*dy editw edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'visible', 'on', ...
  'Tag','resFnString');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~ notes section
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uicontrol('Parent',H0, ...
  'Units','normalized', ...   
  'HorizontalAlignment','left', ...
  'Position',[la3 ba1-16*dy la2-la1+editw edith], ...
  'FontSize',fsz, ...
  'FontWeight','bold', ...
  'BackgroundColor',[0.6 0.6 0.6], ...
  'Style','text', ...
  'String','~~~~~ notes', ...
  'Tag','headline2');

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position', [la3 ba1-18*dy la2-la1+editw 2*edith], ... 
  'FontSize',fsz, ...   
  'BackgroundColor',[1 1 1], ...
  'Style','edit', ...
  'Max',10,...
  'Min',1,...
  'visible', 'on', ...
  'Tag','notesString');


% buttons
uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1 ba1-24*dy butt2x butt2y], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','pushbutton',...  
  'String','load from file', ...
  'TooltipString','',...
  'Tag','loadParBttn', ...
  'callback', {funcH,{'readOptionsFromFile'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+1*butt2x ba1-24*dy butt2x butt2y], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','pushbutton',...  
  'String','save to file', ...
  'TooltipString','',...
  'Tag','saveParBttn', ...
  'callback', {funcH,{'writeOptions2File'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+2*butt2x  ba1-24*dy butt2x butt2y], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','pushbutton',...  
  'String','apply', ...
  'TooltipString','',...
  'Tag','OKParBttn', ...
  'callback', {funcH,{'readOptionsFromGui'}});  

uicontrol('Parent',H0, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[la1+3*butt2x  ba1-24*dy butt2x butt2y], ...
  'FontSize',fsz, ...
  'Fontweight','bold',...
  'style','pushbutton',...  
  'String','apply & close', ...
  'TooltipString','',...
  'Tag','OKPar&CloseBttn', ...
  'callback', {funcH,{'readOptionsFromGui','closeOptionsDialog'}});  
