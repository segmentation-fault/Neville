%   Copyright (C) 2017  Antonio Franco
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

classdef Neville
    %Neville: class based on J. N. Lyness and C. B. Moler, “Van der monde systems and numerical
    %differentiation,” Numerische Mathematik, vol. 8, pp. 458–464, Aug 1966.
    properties
      toDouble;
    end
    methods
        function obj = Neville(varargin)
            obj.toDouble = 1;
            if numel(varargin)
                obj.toDouble = varargin{1};
            end
        end
        
        function v=at(~,V,k)
            assert(k > -1);
            v = V(k+1);
        end
        
        function j = idx(~,k)
            assert(k > -1);
            j = k+1;
        end
        
        function C = update(obj,k,p,x,fxk,C1)
            C = C1;
            xk = obj.at(x,k);
            m = k*(k+3)/2;
            C(obj.idx(m)) = fxk;
            for d=1:k
                xkd = xk - obj.at(x,k-d);
                for s=0:min([d,p])
                    m = m-1;
                    n = m+d;
                    if s==0
                        C(obj.idx(m)) = obj.at(C,n) + xk*(obj.at(C,n-k-1)-obj.at(C,n))/xkd;
                    elseif s==d
                        C(obj.idx(m)) = (obj.at(C,n+1) - obj.at(C,n-k))/xkd;
                    else
                        C(obj.idx(m)) = obj.at(C,n) ...
                            + (xk*(obj.at(C,n-k-1)-obj.at(C,n))+(obj.at(C,n+1)-obj.at(C,n-k)))/xkd;
                    end
                end
                if d>p
                    m = m - d + p;
                end
            end
        end
        
        function [df, r] = dip(obj,p,kmax,f,x,eps)
            %Calculates the df(s) = (s-1)-th derivative of f at O , for s = 1 .. p+1. Uses
            % k-th degree polynomial for some k satisfying p <= k <= kmax. If the rela-
            % tive change in df(s) from k - 1 to k is less than eps then this de-
            % termines k and success is true. Otherwise k = kmax and success is false.
            % p <= length(x) <= kmax
            % df = array of derivatives at 0
            % r = exit flag
            % p = maximum derivative to calculate
            % kmax = maximum degree of the interpolating polynomial
            % f = function
            % x = array of values around 0
            % eps = desired tolerance
            if obj.toDouble
                C = zeros(1,kmax*(kmax+3)/2+1);
                df = zeros(1,p+1);
            else
                C = sym.empty;
                for i=0:kmax*(kmax+3)/2
                    C(obj.idx(i)) = sym(0);
                end
                df = sym.empty;
                for i=0:p
                    df(obj.idx(i)) = sym(0);
                end
            end
            for k=0:kmax
                C = obj.update(k,p,x,f(obj.at(x,k)),C);
                if k < p
                    continue;
                end
                r = 1;
                for s=0:p
                    if r
                        r = abs(obj.at(C,k-s)-obj.at(df,s)) <= eps*abs(obj.at(C,k-s));
                    end
                    df(obj.idx(s)) = obj.at(C,k-s);
                end
                if r
                    break;
                end
            end
            for s=1:p
                df(obj.idx(s)) = factorial(s)*obj.at(df,s);
            end
        end        
    end    
end