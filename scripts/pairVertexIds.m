clear all;

% This matlab code search for opposite id pairs on a cell surface
% if the cell is in the origin, two vertex are oppisite, when the sum 
% of their coordinates almost zero (0,1,-1)+(0,-1,1)= 0 
% how to get the input file: 
%   running a simulation with only one cell
%   load the zero time result into paraview
%   in paraview use only the vertex id pointArray
%   Save the data into a csv file : File/save data
%   Opened csv file header should be looks like this: "Vertex Id,"Points:0","Points:1","Points:2""
%   The centers can be find in tmp/csv/00000000

%% Initialize variables.
filename = 'path/to/saved/csvFilefromParaview.csv';
centerx = 30.0; 
centery = 30.0;
centerz = 30.0;
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Create output variable
wbc = table(dataArray{1:end-1}, 'VariableNames', {'VertexId','Points0','Points1','Points2'});

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
wbcN = wbc;
% substract the center, so we have the WBC in the origin
wbcN.Points0 = wbc.Points0 - centerx;
wbcN.Points1 = wbc.Points1 - centery;
wbcN.Points2 = wbc.Points2 - centerz;
pairs = zeros(2,length(wbcN.Points1)/2);
k = 1;
error = 5e-3;
for i = 1:length(wbcN.Points0)
    for j = i:length(wbcN.Points0)
        if abs(wbcN.Points0(i)+ wbcN.Points0(j)) < error 
            if abs(wbcN.Points1(i)+ wbcN.Points1(j)) < error
                if abs(wbcN.Points2(i) + wbcN.Points2(j)) < error
                    pairs(1,k) = wbcN.VertexId(i);
                    pairs(2,k) = wbcN.VertexId(j);
                    k = k+1;
                end
            end
            
        end
    end
end
missed = [];
for m=0:length(wbcN.Points0)-1
    contains = ismember(m,pairs);
    if contains == 0
        missed = [missed, m];
    end
end
if isempty(missed) == 0
    fprintf("Some point somehow skipped. Check the maximum error\n");
    return;
end

formatSpec = '<Edge> %d %d </Edge>\n';
fileID = fopen('edges.txt','w');
fprintf(fileID, formatSpec,pairs);
fclose(fileID);