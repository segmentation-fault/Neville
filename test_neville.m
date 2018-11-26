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

%Test for the Neville algorithm against the symbolic differentiation in
%Matlab

clc
clear
close all
reset(symengine)

tic

syms t real;

nPoints = 100; %100 points around 0
p = 6; %Calculate all the derivatives until the 6th
x = -1; %Deriving in -1

%Defining the function
F = 6*exp(-6*t)+t+log(0.2*t+4)-2^t-erf(t);
f = matlabFunction(F);

%Calculating the derivative numerically
N = Neville();
N.toDouble = false; %Using VPA
[df, succ] = N.deriveAt(p,f,x,1e-6,linspace(-0.1,0.1,nPoints));
N.toDouble = true; %Not using VPA
[df1, succ1] = N.deriveAt(p,f,x,1e-6,linspace(-0.1,0.1,nPoints));

%Symbolical derivation at x
Fv = sym.empty(0,p+1);
Fv(1) = subs(F,t,0);
for i=1:p
    Fv(i+1) = subs(diff(F,t,i),t,x);
end

%Plot
figure;
plot(0:p,df,'rx:',0:p,df1,'k*-.',0:p,Fv,'b+--');
legend({'Neville (VPA)','Neville (no VPA)','Symbolic'});
xlabel('order of derivative')
title(['Derivatives in t = ' num2str(x) ' of $$f(t) = ' latex(F) '$$'],'interpreter','latex')

elapsedTime = toc;
s = seconds(elapsedTime);
s.Format = 'hh:mm:ss.SSS';
disp('Elapsed time (hms): ');
disp(s);