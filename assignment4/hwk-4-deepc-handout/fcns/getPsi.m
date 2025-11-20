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

function [psi] = getPsi(A,B,Tf,type)


switch type
    case {'LTI', 'lti'}

        [n,m] = size(B);
        psi = zeros(n*Tf, m*Tf);
        
        %first column
        psicol = zeros(n*Tf, m);
        psicol(n+1:2*n,:) = B;
        
        for ii = 2:Tf-1
            psicol(ii*n+1:(ii+1)*n,:)  = A*psicol((ii-1)*n+1:ii*n,:);
        end
        
        % there are tf columns, and the last is all zeros (so skip it). we
        % have generated the first one in psicol
        
        for jj = 2:Tf-1
            psi((jj-1)*n+1:end,(jj-1)*m+1:jj*m) = psicol(1:n*Tf-(jj-1)*n,:);
        end
        
        % actually insert the first column lol
        
        psi(:,1:m) = psicol;

    case {'LTV', 'ltv'}
        
        [n,m] = size(B(:,:,1));
        psi = zeros(n*Tf, m*Tf);
        
        %LTV case is done row-by-row by computing the upsilon matrices        
        for ii = 2:Tf
            upsilon = getUpsilon(A,B,ii);
            [~, ups_col] = size(upsilon);
            psi((ii-1)*n+1:ii*n,1:ups_col) = upsilon;
        end

end

