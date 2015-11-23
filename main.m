%% DEMO
% Results of the relaxation labelling algorithm on three different matrices
% of probabilities and their compatibility coefficients
% author: Loubna Rizqi
% date: 21/11/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Demo 1
% Probability matrix: pixels being edges
% P = [ 0   0   0   0 0
%       0 0.1 0.1 0.1 0
%       0 0.1 0.9   0 0
%       0 0.1   0   0 0 
%       0   0   0   0 0 ]
%%
init_1 = zeros(5,5);
init_1(2:4,2) = 0.1;
init_1(2,2:4) = 0.1;
init_1(3,3) = 0.9;

% Compatibility coefficients
% column & row 1: edge
% column & row 2: non edge
r_1 = ones(2,2);

%% Demo 2
% Probability matrix: pixels being edges
% P = [ 0 0 0 0 0
%       0 0 0 0 0
%       0 0 1 0 0
%       0 0 0 0 0 
%       0 0 0 0 0 ]
%%
init_2 = zeros(5,5);
init_2(3,3)=1;

% Compatibility coefficients
% column & row 1: edge
% column & row 2: non edge
r_2 = ones(2,2);
r_2(1,1)=2;

%% Demo 3
% Probability matrix: pixels being edges
% P = [ 0   0   0   0 0
%       0 0.1 0.1 0.1 0
%       0 0.1   1 0.1 0
%       0 0.1 0.1 0.1 0 
%       0   0   0   0 0 ]
%%
init_3 = zeros(5,5);
init_3(2:4,2:4)=0.1;
init_3(3,3)=1;

%% Number of iterations
n=2; % Let's start with with 2

%% Results
disp('--------------------------------------------')
disp('Relaxation labelling results for demo 1:')
disp('--------------------------------------------')
display(relax_labelling(init_1,r_1,n))
disp('--------------------------------------------')
disp('Relaxation labelling results for demo 2:')
disp('--------------------------------------------')
display(relax_labelling(init_2,r_2,n))
disp('--------------------------------------------')
disp('Relaxation labelling results for demo 3:')
disp('--------------------------------------------')
display(relax_labelling(init_3,r_2,n))

