function output = relax_labelling(Pe, r, n)
%% Applies the relaxation labelling algorithm to an initial
% 5x5 matrix of probabilities and compatibility coefficients
% INPUT
% Pe: Initial probabilities of each pixel being an edge
%     (borders are zeros)
% r: Compatibility coefficients
% n: number of iterations
% OUTPUT
% new matrix of probabilities of each pixel being an edge
%
% author: Loubna Rizqi
% date: 21/11/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Probability matrix: pixels not being edges
Pne = ones(5,5) - Pe;

[x,y] = find(Pe~=0); % Probabilities equal to zero remain unchanged
% Building the Weighting Coefficients
C = zeros(5,5,5,5);
for k=1:length(x)
    C(x(k),y(k),x(k)-1:x(k)+1,y(k)-1:y(k)+1)=1;
    C(x(k),y(k),x(k),y(k))=0;
end

% Q matrix, the 2 corresponds to e and ne
% Q(:,:,1): label is "edge"
% Q(:,:,2): label is "non edge"
Q = zeros(5,5,2);
% Applying the relaxation schema
for k=1:n
    for i=2:4
        for j=2:4
            Q(i,j,1) = C(i,j,i-1,j)*(r(1,1)*Pe(i-1,j)+r(1,2)*Pne(i-1,j))+...
                C(i,j,i,j+1)*(r(1,1)*Pe(i,j+1)+r(1,2)*Pne(i,j+1))+...
                C(i,j,i-1,j-1)*(r(1,1)*Pe(i-1,j-1)+r(1,2)*Pne(i-1,j-1))+...
                C(i,j,i-1,j+1)*(r(1,1)*Pe(i-1,j+1)+r(1,2)*Pne(i-1,j+1))+...
                C(i,j,i,j-1)*(r(1,1)*Pe(i,j-1)+r(1,2)*Pne(i,j-1))+...
                C(i,j,i+1,j)*(r(1,1)*Pe(i+1,j)+r(1,2)*Pne(i+1,j))+...
                C(i,j,i+1,j+1)*(r(1,1)*Pe(i+1,j+1)+r(1,2)*Pne(i+1,j+1))+...
                C(i,j,i+1,j-1)*(r(1,1)*Pe(i+1,j-1)+r(1,2)*Pne(i+1,j-1));
            Q(i,j,2) = C(i,j,i-1,j)*(r(2,1)*Pe(i-1,j)+r(2,2)*Pne(i-1,j))+...
                C(i,j,i,j+1)*(r(2,1)*Pe(i,j+1)+r(2,2)*Pne(i,j+1))+...
                C(i,j,i-1,j-1)*(r(2,1)*Pe(i-1,j-1)+r(2,2)*Pne(i-1,j-1))+...
                C(i,j,i-1,j+1)*(r(2,1)*Pe(i-1,j+1)+r(2,2)*Pne(i-1,j+1))+...
                C(i,j,i,j-1)*(r(2,1)*Pe(i,j-1)+r(2,2)*Pne(i,j-1))+...
                C(i,j,i+1,j)*(r(2,1)*Pe(i+1,j)+r(2,2)*Pne(i+1,j))+...
                C(i,j,i+1,j+1)*(r(2,1)*Pe(i+1,j+1)+r(2,2)*Pne(i+1,j+1))+...
                C(i,j,i+1,j-1)*(r(2,1)*Pe(i+1,j-1)+r(2,2)*Pne(i+1,j-1));
        end
    end
    
    for i=1:5
        for j=1:5
            if Pe(i,j)~= 0
                % Calculating the new matrices
                sum = Pe(i,j)*Q(i,j,1)+Pne(i,j)*Q(i,j,2);
                Pe(i,j) = Pe(i,j)*Q(i,j,1)/sum;
            end
        end
    end
    Pne = ones(5,5)-Pe;
end
            
output = Pe;

end
 
