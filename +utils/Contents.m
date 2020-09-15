% +UTILS  Package for utilities in Kalman-Thesis repository 
% MATLAB Version 9.7 (R2019b Update 5) 29-Aug-2020
%
% Functions
%  addAreaToAxes                 - Add label to axes indicating 'RFA' or 'S1' in a principled way
%  addConnectingLine             - Add (annotated) line connecting point pairs <X,Y>
%  addLabelsToAxes               - Add default labels to axes
%  addLegendToAxes               - Add legend to axes depending on fields of params
%  addPeakLabels                 - Add labels denoting <x,y> value at list of peaks
%  addSolenoidToAxes             - Add indicator of solenoid strike to axes
%  addStimInfoToAxes             - Add information about stimulus-type, timing, and if ICMS occurred on this channel.
%  addTextToAxes           		 - Add text to a specified location on an axes
%  addTimeIndicatorToAxes  		 - Add graphics indicator for some event to an axes
%  addTypeToAxes           		 - Add label to axes indicating 'Solenoid (only)', 'ICMS (only)', or 'Solenoid+ICMS' in a principled way
%  checkXYLabels                 - Check X-Y axes for correct labels at end of label strings
%  expAI                         - Export figure in appropriate format for Adobe Illustrator
%  formatDefaultAxes             - Apply default axes settings preferred by MM
%  formatDefaultFigure           - Apply default figure settings preferred by MM
%  formatDefaultLabel            - Apply default font settings preferred by MM
%  getFigAx                      - Return figure and axes handles given parameters struct
%  parseTitle                    - Parses title from filter input arguments
%  plotCoherence                 - Helper function to plot coherence (from example R2017a)
%  plotSpikeRaster               - Create raster plot from binary spike data or spike times
%
% <strong>Functions</strong>
%  fastsmooth                    - Smooths vector X
%  getCB95                       - Return 95% confidence bounds
%  getKinematicStates            - Convert initial kinematic table P into states for cleaning
%  getOpt                        - Process paired optional arguments as `prop`,val1 
%  getPathTo                     - Return output from uigetdir basically
%  getPredictionMask             - Return mask for matched samples for fitting prediction
%  makeKey                       - Utility to make random alphanumeric key-string for naming a "row"
%  parseNamedVariable            - Parse optional inputs and update a named variable
%  sec2string                    - Take seconds (double) and return time string for hours, minutes, and seconds