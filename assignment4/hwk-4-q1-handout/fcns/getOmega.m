% MIT License
% 
% Copyright (c) 2021, ETH Zurich, [Mathias Hudoba de Badyn]
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
%     Running this software requires a local Gurobi installation. Academic
%     licenses, at the time of writing, are available at no cost at
%     http://www.gurobi.com

function [omega] = getOmega(A,Tf,type)
% outputs the n*Tf by n matrix Omega, whose n by n blocks are the state
% transition matrices, with time increasing as the blocks go down


switch type
    case {'LTI', 'lti'}

        [n,~] = size(A);
        omega = [eye(n,n); zeros(n*(Tf-1), n)];
        
        % A is static in the LTI case, just multiply each block by A
        for ii = 2:Tf
            omega((ii-1)*n+1:ii*n,:) = A*omega((ii-2)*n+1:(ii-1)*n,:);
        end

    case {'LTV', 'ltv'}
        
        [n,~] = size(A(:,:,1));
        omega = [eye(n,n); zeros(n*(Tf-1), n)];
        
        % put the state transition matrix in for each block of omega
        for ii = 2:Tf
            omega((ii-1)*n+1:ii*n,:) = A(:,:,ii-1)*omega((ii-2)*n+1:(ii-1)*n,:);
        end

end
