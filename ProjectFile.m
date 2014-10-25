%**************************************************************************
function ProjectFile(action)
%**************************************************************************
% 
% File name: ProjectFile.m
% Author: Khurom H. Kiyani [k.kiyani@warwick.ac.uk]
% File created: 13th March 2014
% Last updated: 21st May 2014
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This function is intended for project organisation. Be sure to be in your
% main project directory (MPD) when you call this function. 
% 
% Input arguments
% 
% 'open':   When launched with this argument, all the functions and scripts 
%           contained within the 'WorkFiles.mat' file in your MPD [which
%           was saved earlier through lauching the 'close' argument in this 
%           function (see next) or through intialising (see notes)] and 
%           are located in the 'files' variable, are opened up in the
%           Matlab Editor. Also the local functions and scripts in the 
%           'pathstuff' variable will be added to the Matlab path.
%   
% 'close':  When launched with this argument, all the scripts and functions 
%           which are open in the Matlab Editor will be saved to a
%           'filenames' variable in 'WorkFiles.mat' in the MPD, ready to be
%           opened again later using the 'open' argument. Also the paths 
%           which contain local functions and scripts and were added 
%           earlier will be removed from the Matlab path. 
% 
% NOTES [Initialisation]:
% -----------------------
% This function expects that a file called 'WorkFiles.mat' is in the MPD,
% and that it contains two variables: 'files' and 'pathstuff'. When you
% create a project, make sure to include these variables in this file
% 'WorkFiles.mat' which is saved to your MPD even if the 'files' variable
% contains nothing initially. 'pathstuff' should be made by going into the
% folder in your MPD which contains your functions and scripts (even if it
% contains nothing initially), and typing in the command:
% 
%   pathstuff=genpath(pwd)
% 
% Then an empty cell variable 'files' should be created by the command:
% 
%   files=cell(1,1)
% 
% Then both of these should be saved to your MPD with the command:
% 
%   save WorkFiles.mat files pathstuff
% 
% Through the course of your project these variables will be populated by 
% the files you open up and the local functions that you will create. 
% 
%**************************************************************************

switch action
   
    % When you first start the project
    
    case 'open'
        
        filesObj=matlab.desktop.editor.getAll;
        files={filesObj.Filename}';
        
        if isempty(files)==0
            
            a=input(['Are you sure? You already have stuff open;', ...
           ' this stuff will be added to it (1=yes, anything else=no): ']);
       
        else
            
            a=1;
            
        end
            
        clear files
        clear filesObj
        
        if a==1
            
            load WorkFiles.mat;
            
            % Add functions & scripts local to this project to Matlab path
            addpath(pathstuff);
        
            % Open up files in the Matlab editor
            for k=1:numel(files)
                open(files{k});
            end
        
            clear k files pathstuff
        
        else
            
            disp('Phew, that was close; you nearly screwed up!')
            
        end
        
        clear a
        
    % When you close the project    
        
    case 'close'
        
        filesObj=matlab.desktop.editor.getAll;
        files={filesObj.Filename}';
        
        if isempty(files)==1
            
            a=input(['Are you sure? You have nothing open;', ...
           ' existing WorkFiles.mat will be replaced with nothing ', ...
           '(1=yes, anything else=no): ']);
       
        else
            
            a=1;
            
        end
            
        clear files
        clear filesObj
        
        if a==1
        
            load WorkFiles.mat;
            clear files
        
            % Remove functions & scripts local to this project from Matlab path
            rmpath(pathstuff);
        
            % Close all open files in the Matlab editor and save them
            filesObj=matlab.desktop.editor.getAll;
            files={filesObj.Filename}';
            save WorkFiles.mat files pathstuff
    
            clear files pathstuff
    
            filesObj.close
            clear filesObj
            
        else
            
            disp('Phew, that was close; you nearly screwed up!')
            
        end
        
        clear a
        
    otherwise
        
        disp('!!!Not a valid argument!!!')
    
end

end