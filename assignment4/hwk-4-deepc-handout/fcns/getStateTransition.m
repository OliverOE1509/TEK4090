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

function [psi] = getStateTransition(A,l,m)

% outputs the state transition matrix psi(l,m) = A_{l-1} * ... * A_{m}
% and with psi(1,1) = eye. Assumes that A has time in third dimension.

if l<m
    warning('l is less than m, i.e. the arguments of the state transition matrix are backwards')
end

[~,~,c] = size(A);

if c>1
    % system is LTV

    psi = eye(size(A(:,:,1)));

    for ii = m:l-1
        psi = A(:,:,ii)*psi;
    end

else
    % system is LTI
    
    psi = eye(size(A(:,:,1)));

    for ii = m:l-1
        psi = A(:,:)*psi;
    end
    
end




end
