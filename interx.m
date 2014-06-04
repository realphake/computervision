% Copyright (c) 2009, Hrishi Shah
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate intersection of three spheres
% Returns NAN if no intersection or some bad condition encountered
% X1, X2, X3 are vectors of centers of three spheres
% r1, r2, r3 are radii of the three spheres
% pos =0 for lower point and 1 for higher point
% Computation Code calculated in (and copied from) Maple
% For questions/suggestions, contact hrishi.shah2002@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=interx(X1,X2,X3,r1,r2,r3,pos)
if(nargin<7), pos=1; end % default value
x1=X1(1); y1=X1(2); z1=X1(3);
x2=X2(1); y2=X2(2); z2=X2(3);
x3=X3(1); y3=X3(2); z3=X3(3);
% convert in coord sys at [x1 y1 z1] oriented same as global
% x2=1; y2=1; z2=0.2; x3=2; y3=1; z3=1; r1=2.5; r2=2.6; r3=2.7; % TEST
x2=x2-x1; y2=y2-y1; z2=z2-z1;
x3=x3-x1; y3=y3-y1; z3=z3-z1;
a=(16*y2^2*z3*y3^2*z2*x3*r1^2*x2-4*y2^3*z3*y3*z2*x3*r1^2*x2+4*y2^3*z3*y3*z2*x3*x2*r3^2-4*y2*y3^3*z2*x2*z3*r1^2*x3+4*y2*y3^3*z2*x2*z3*r2^2*x3+16*z2*x3^2*x2^2*r1^2*y3*y2*z3-4*z2*x3^3*x2*r1^2*y3*y2*z3+4*z2*x3^3*x2*y2*z3*r2^2*y3-4*x2^3*z3*x3*y2*r1^2*z2*y3+4*x2^3*z3*x3*y2*r3^2*z2*y3-4*y2*z3*z2^3*y3*x3*r1^2*x2+4*y2*z3*z2^3*y3*x3*x2*r3^2+8*y2*z3^2*z2^2*y3*x2*r1^2*x3-4*y2*z3^2*z2^2*y3*x2*r2^2*x3+4*x2^2*y3^2*z2*y2^2*z3*x3^2-4*x2^4*z3^2*x3^2*y3*y2-2*z2^2*x3^4*x2^2*y2^2+2*z2^2*x3^3*y2^2*r1^2*x2-2*z2^2*x3^3*y2^2*x2*r3^2+2*z2^2*x3^5*x2*y2*y3+2*x2^5*z3^2*x3*y3*y2+2*z2^3*y3^2*y2^2*z3*x3^2+2*z2^3*y3^2*x2^2*z3*x3^2+2*z2^3*y3^2*x2^2*z3*r1^2-2*z2^3*y3^2*x2^2*z3*r3^2+2*x2^2*y3^4*z2*y2^2*z3-4*x2^2*y3^2*z2^2*x3^2*y2^2+2*x2^4*y3^2*z2*z3*x3^2-2*x2^2*y3^3*z2^2*y2*x3^2+2*x2^2*y3^3*z2^2*y2*z3^2+2*x2^3*y3^2*z2^2*x3*z3^2+2*x2^2*y3^2*z2*y2^2*z3^3+2*x2^2*y3^3*z2^2*y2*r1^2+2*x2^2*y3^2*z2^2*x3^2*r2^2+2*x2^4*y3^2*z2*z3*r1^2-2*x2^4*y3^2*z2*z3*r3^2-2*x2^2*y3^3*z2^2*y2*r3^2+2*x2^3*y3^2*z2^2*x3*r1^2-2*x2^3*y3^2*z2^2*x3*r3^2+2*y2^4*z3*x3^2*y3^2*z2+2*y2^2*z3*x3^4*z2*x2^2-4*y2^2*z3^2*x3^2*x2^2*y3^2+2*y2^3*z3^2*x3^2*z2^2*y3+2*y2^2*z3^2*x3^3*x2*z2^2-2*y2^3*z3^2*x3^2*x2^2*y3+2*y2^3*z3^2*x3^2*r1^2*y3+2*y2^2*z3*x3^4*z2*r1^2-2*y2^2*z3*x3^4*z2*r2^2+2*y2^2*z3^2*x3^2*x2^2*r3^2-2*y2^3*z3^2*x3^2*r2^2*y3+2*y2^2*z3^2*x3^3*x2*r1^2-2*y2^2*z3^2*x3^3*x2*r2^2-2*y2^2*z3^2*y3^2*x2^3*x3-4*y2^4*z3^2*y3^2*x2*x3+2*y2^2*z3^2*y3^2*x2^2*r3^2+4*y2^3*z3^2*y3*x2^3*x3+2*y2^5*z3^2*y3*x2*x3+4*y2*y3^3*z2^2*x3^3*x2+2*y2*y3^5*z2^2*x3*x2-2*y2^2*y3^2*z2^2*x3^3*x2-4*y2^2*y3^4*z2^2*x3*x2+2*y2^2*y3^2*z2^2*x3^2*r2^2-4*z2^2*x3^4*x2^2*y2*y3+2*z2*x3^2*x2^2*y2^2*z3^3+2*z2*x3^2*y2^4*z3*r1^2+2*z2^2*x3^2*y2^3*r1^2*y3-2*z2*x3^2*y2^4*z3*r3^2-2*z2^2*x3^2*y2^3*r3^2*y3-z2^4*y3^2*x3^2*x2^2-z2^4*y3^2*x3^2*y2^2+2*z2^3*y3^4*x2^2*z3+2*z2^3*y3^2*x2^2*z3^3+2*x2^2*y3^5*z2^2*y2-2*x2^2*y3^4*z2^2*y2^2-2*x2^4*y3^2*z2^2*x3^2+2*x2^3*y3^2*z2^2*x3^3+2*x2^4*y3^4*z2*z3+2*x2^3*y3^4*z2^2*x3+2*x2^2*y3^4*z2^2*r2^2-2*y2^4*z3^2*x3^2*y3^2+2*y2^5*z3^2*x3^2*y3+2*y2^4*z3*x3^4*z2+2*y2^2*z3^2*x3^3*x2^3-2*y2^2*z3^2*x3^4*x2^2+2*y2^4*z3^2*x3^3*x2+2*y2^2*z3*x3^4*z2^3-y2^2*z3^4*x3^2*x2^2+2*y2^4*z3^2*x3^2*r3^2-2*y2^2*z3^2*y3^4*x2^2+2*y2^3*z3^2*y3^3*x2^2-y2^2*z3^4*y3^2*x2^2-y2^4*z3^2*y3^2*x2^2+2*y2^3*y3^3*z2^2*x3^2-y2^2*y3^4*z2^2*x3^2-2*y2^4*y3^2*z2^2*x3^2-z2^4*y3^4*x2^2-2*x2^4*y3^4*z2^2-2*y2^4*z3^2*x3^4-y2^4*z3^4*x3^2-2*z2^2*x3^4*y2^4-z2^4*x3^4*y2^2-2*x2^4*z3^2*y3^4-x2^4*z3^4*y3^2-z3^2*y2^6*x3^2-x3^6*z2^2*y2^2-z3^2*x2^6*y3^2-2*y2^4*x3^4*y3^2+2*y2^4*x3^4*r3^2+2*y2^5*x3^4*y3-y2^4*x3^2*r1^4-y2^4*x3^2*y3^4+2*y2^5*x3^2*y3^3-y2^4*x3^2*r3^4-y2^6*x3^2*y3^2-y2^2*x3^4*r1^4-y2^2*x3^4*x2^4+2*y2^2*x3^5*x2^3-y2^2*x3^4*r2^4-y2^2*x3^6*x2^2-2*y2^4*x3^4*x2^2+2*y2^4*x3^4*r2^2+2*y2^4*x3^5*x2+2*x2^4*y3^5*y2-2*x2^4*y3^4*y2^2+2*x2^4*y3^4*r2^2-x2^2*y3^4*r1^4-x2^2*y3^6*y2^2+2*x2^2*y3^5*y2^3-x2^2*y3^4*y2^4-x2^2*y3^4*r2^4-x2^4*y3^2*r1^4-x2^6*y3^2*x3^2+2*x2^5*y3^2*x3^3-x2^4*y3^2*x3^4-x2^4*y3^2*r3^4+2*x2^5*y3^4*x3-2*x2^4*y3^4*x3^2+2*x2^4*y3^4*r3^2-y3^6*z2^2*x2^2-y2^4*x3^6-y2^6*x3^4-x2^6*y3^4-x2^4*y3^6-2*z2*x3^2*r1^2*y2^2*z3*r3^2-2*z2*y3^2*r2^2*x2^2*z3*r1^2+2*z2*y3^2*r2^2*x2^2*z3*r3^2+2*x2^4*y3^2*z2*z3^3+2*y2*r1^4*z2^2*y3*x3*x2+2*x2^2*y3^2*z2*y2^2*z3*r1^2-8*x2^2*y3^3*z2*r1^2*y2*z3-2*x2^2*y3^2*z2*y2^2*z3*r3^2-8*x2^3*y3^2*z2*z3*r1^2*x3+2*y2^2*z3*x3^2*r1^2*y3^2*z2-8*y2^3*z3*x3^2*r1^2*z2*y3-2*y2^2*z3*x3^2*z2*y3^2*r2^2-8*y2^2*z3*x3^3*z2*r1^2*x2-4*y2^2*z3^2*y3^2*x2*z2^2*x3-4*y2^2*z3^2*y3^2*x2*r1^2*x3+4*y2^2*z3^2*y3^2*x2*r2^2*x3-4*y2^3*z3*y3*z2*x3^3*x2-4*y2^3*z3*y3^3*z2*x3*x2-4*y2^3*z3^3*y3*z2*x3*x2+4*y2^3*z3^2*y3*x2*z2^2*x3-4*y2^3*z3^2*y3*x2*r2^2*x3-4*y2*y3^3*z2*x2^3*z3*x3+4*y2*y3^3*z2^2*x3*x2*z3^2-4*y2*y3^3*z2^3*x2*z3*x3-4*y2*y3^3*z2^2*x3*x2*r3^2-4*y2^2*y3^2*z2^2*x3*r1^2*x2+4*y2^2*y3^2*z2^2*x3*x2*r3^2-4*z2^2*x3^2*x2^2*y2*z3^2*y3+2*z2*x3^2*x2^2*y2^2*z3*r1^2-4*z2^2*x3^2*x2^2*y2*r1^2*y3-2*z2*x3^2*x2^2*y2^2*z3*r3^2+4*z2^2*x3^2*x2^2*y2*r3^2*y3-4*z2^3*x3^3*x2*y2*z3*y3+4*z2^2*x3^3*x2*y2*z3^2*y3-4*z2*x3^3*x2^3*y3*y2*z3-4*z2^2*x3^3*x2*y2*r3^2*y3+4*x2^3*z3^2*x3*y2*z2^2*y3-4*x2^3*z3^3*x3*y2*z2*y3-4*x2^3*z3^2*x3*y2*r2^2*y3+2*x2^2*z3*x3^2*r1^2*y3^2*z2-4*x2^2*z3^2*x3^2*r1^2*y3*y2-2*x2^2*z3*x3^2*z2*y3^2*r2^2+4*x2^2*z3^2*x3^2*y2*r2^2*y3-4*y2*z3^3*z2^3*y3*x3*x2+2*y2*z3^2*z2^4*y3*x2*x3+2*y2*z3^4*z2^2*y3*x3*x2-2*r1^2*y3^2*z2*x2^2*z3*r3^2-2*y2^2*z3*r1^2*z2*x3^2*r2^2+2*r1^4*y3*y2*z3^2*x2*x3+2*z2^2*x3^5*y2^2*x2+2*z2^2*x3^4*y2^3*y3+2*z2*x3^2*y2^4*z3^3+2*z2^2*x3^4*y2^2*r2^2-z2^2*x3^4*x2^2*y3^2+2*x2^5*z3^2*x3*y3^2-x2^4*z3^2*x3^2*y2^2-2*x2^4*z3^2*x3^2*y3^2+2*x2^4*z3^2*y3^3*y2+2*x2^4*z3^2*y3^2*r3^2-2*y2^2*x3^4*z2^2*y3^2-2*z2^2*x3^2*x2^2*y3^4-2*x2^2*z3^2*y2^4*x3^2-2*x2^4*y3^2*y2^2*z3^2+2*y2^2*z3^3*z2^3*x3^2-z3^2*y2^2*r1^4*x3^2-z3^2*y2^2*z2^4*x3^2-z3^2*y2^2*r2^4*x3^2-2*z3^2*y2^4*x3^2*z2^2+2*z3^2*y2^4*x3^2*r2^2-2*x3^4*z2^2*y2^2*z3^2+2*x3^4*z2^2*y2^2*r3^2-x3^2*z2^2*y2^2*r1^4-x3^2*z2^2*y2^2*z3^4-x3^2*z2^2*y2^2*r3^4-2*z3^2*x2^4*y3^2*z2^2+2*z3^2*x2^4*y3^2*r2^2-z3^2*x2^2*r1^4*y3^2-z3^2*x2^2*z2^4*y3^2-z3^2*x2^2*r2^4*y3^2+2*y2^3*x3^4*r1^2*y3-2*y2^3*x3^4*x2^2*y3-2*y2^3*x3^4*r2^2*y3+2*y2^4*x3^3*r1^2*x2-2*y2^4*x3^3*x2*y3^2-2*z2^2*x3^2*x2^2*y3^2*z3^2+2*z2^2*x3^2*x2^2*y3^2*r3^2-2*x2^2*z3^2*y2^2*x3^2*z2^2+2*x2^2*z3^2*y2^2*x3^2*r2^2+2*x2^2*y3^2*y2^2*z3^2*r2^2+2*y2^2*z3^3*z2*x3^2*r1^2-2*y2^2*z3^3*z2*x3^2*r2^2+2*z2^3*x3^2*y2^2*z3*r1^2-2*z2^3*x3^2*y2^2*z3*r3^2+2*x2^2*z3^3*r1^2*y3^2*z2-2*x2^2*z3^3*z2*y3^2*r2^2+2*r1^4*y3^2*z2*x2^2*z3+2*y2^2*z3*r1^4*z2*x3^2+2*x2^2*z3*y3^4*r1^2*z2+2*x2^2*z3^2*y3^3*r1^2*y2-2*x2^2*z3*y3^4*z2*r2^2-2*x2^2*z3^2*y3^3*y2*r2^2+2*x2^3*z3^2*y3^2*r1^2*x3-2*x2^3*z3^2*y3^2*r2^2*x3-2*y2^2*z3^2*z2^2*y3^2*x2^2-2*y2^2*x3^2*z2^2*y3^2*z3^2+2*y2^2*x3^2*z2^2*y3^2*r3^2-4*y2*z3^2*z2^2*y3*x3*x2*r3^2-4*y2*z3^3*z2*y3*x2*r1^2*x3+4*y2*z3^3*z2*y3*x2*r2^2*x3-4*r1^4*y3*y2*z3*z2*x3*x2-4*r1^2*y3*y2*z3^2*x2*r2^2*x3-4*y2*r1^2*z2^2*y3*x3*x2*r3^2+4*r1^2*y3*y2*z3*z2*x3*x2*r3^2+4*y2*r1^2*z2*y3*x2*z3*r2^2*x3+2*z2*x3^2*r2^2*y2^2*z3*r3^2+2*y2*z3^2*r2^4*y3*x2*x3+2*y2*r3^4*z2^2*y3*x3*x2+4*y2^2*x3*x2*y3^2*r1^2*r3^2+4*y2^2*x3*x2*y3^2*r1^2*r2^2-4*y2^2*x3*x2*y3^2*r3^2*r2^2+4*y2*x3^2*x2^2*y3*r1^2*r2^2+4*y2*x3^2*x2^2*y3*r1^2*r3^2-4*y2*x3^2*x2^2*y3*r2^2*r3^2-4*y2^3*x3*x2*y3*z3^2*r3^2-4*y2*x3*x2*y3^3*r1^2*r2^2-4*y2^3*x3*x2*y3*r1^2*r3^2-4*y2*x3*x2*y3^3*z2^2*r2^2-4*y2*x3*x2^3*y3*r1^2*r3^2-4*y2*x3^3*x2*y3*r1^2*r2^2-4*y2*x3^3*x2*y3*z2^2*r2^2-4*y2*x3*x2^3*y3*z3^2*r3^2-4*z3^2*y2^2*r1^2*x3^2*z2^2+2*z3^2*y2^2*r1^2*x3^2*r2^2+2*z3^2*y2^2*z2^2*x3^2*r2^2+2*x3^2*z2^2*y2^2*z3^2*r3^2+2*x3^2*z2^2*y2^2*r1^2*r3^2-4*z3^2*x2^2*r1^2*y3^2*z2^2+2*z3^2*x2^2*r1^2*y3^2*r2^2+2*z3^2*x2^2*z2^2*y3^2*r2^2-2*y2^3*x3^2*r1^2*y3*r3^2-2*y2^3*x3^2*x2^2*y3*r1^2+2*y2^3*x3^2*x2^2*y3*r3^2-2*y2^3*x3^2*r1^2*r2^2*y3+2*y2^3*x3^2*r3^2*r2^2*y3-2*y2^2*x3^3*r1^2*x2*r2^2-2*y2^2*x3^3*r1^2*x2*y3^2-2*y2^2*x3^3*r1^2*x2*r3^2+2*y2^2*x3^3*r2^2*x2*y3^2+2*y2^2*x3^3*r2^2*x2*r3^2+2*y2^2*x3^2*r1^2*y3^2*r2^2+4*y2^2*x3^2*x2^2*y3^2*r2^2+2*y2^2*x3^2*r1^2*x2^2*r3^2+4*y2^2*x3^2*x2^2*y3^2*r3^2-8*y2^3*x3^3*r1^2*x2*y3-2*x2^2*y3^3*r1^2*y2*x3^2-2*x2^2*y3^3*r1^2*y2*r3^2-2*x2^2*y3^3*y2*r1^2*r2^2+2*x2^2*y3^3*y2*x3^2*r2^2+2*x2^2*y3^3*y2*r3^2*r2^2-2*x2^3*y3^2*r1^2*y2^2*x3-2*x2^3*y3^2*r1^2*r2^2*x3-2*x2^3*y3^2*r1^2*x3*r3^2+2*x2^3*y3^2*y2^2*x3*r3^2+2*x2^3*y3^2*r2^2*x3*r3^2+2*x2^2*y3^2*y2^2*r1^2*r3^2-4*y2*z3*r2^2*y3*z2*x3*x2*r3^2-4*x2^4*y3^2*r1^2*x3^2+2*x2^4*y3^2*r1^2*r3^2+2*x2^3*y3^2*r1^2*x3^3+2*x2^4*y3^2*x3^2*r2^2-2*x2^5*y3^2*x3*r3^2-2*x2^3*y3^2*r2^2*x3^3+2*x2^4*y3^2*x3^2*r3^2-y3^2*z2^2*r1^4*x2^2-y3^2*z2^2*x2^2*z3^4-y3^2*z2^2*x2^2*r3^4-2*y3^4*z2^2*x2^2*z3^2+2*y3^4*z2^2*x2^2*r3^2+4*y2^3*x3*x2^3*y3^3+4*y2^3*x3^3*x2*y3^3+2*y2*x3*x2^5*y3^3+2*y2^3*x3^5*x2*y3+2*y2^3*x3*x2*y3^5-4*y2^4*x3*x2*y3^4+2*y2^5*x3*x2*y3^3+2*y2*x3^3*x2^5*y3-4*y2*x3^4*x2^4*y3+2*y2^5*x3^3*x2*y3+2*y2*x3^5*x2^3*y3+2*y2*x3*x2^3*y3^5+4*y2^3*x3^3*x2^3*y3+4*y2*x3^3*x2^3*y3^3-2*y2^4*x3^3*x2*r3^2-2*y2^5*x3^2*r3^2*y3+2*y2^4*x3^2*y3^2*r2^2+2*y2^3*x3^2*r1^2*y3^3+2*y2^3*x3^2*r1^4*y3-4*y2^4*x3^2*r1^2*y3^2-3*y2^4*x3^2*x2^2*y3^2+2*y2^4*x3^2*r1^2*r3^2+2*y2^5*x3^2*r1^2*y3+2*y2^4*x3^2*y3^2*r3^2-2*y2^3*x3^2*y3^3*r2^2-y2^2*x3^2*r1^4*y3^2-3*y2^2*x3^2*x2^4*y3^2-y2^2*x3^2*r2^4*y3^2-y2^2*x3^2*r1^4*x2^2-3*y2^2*x3^2*x2^2*y3^4-y2^2*x3^2*x2^2*r3^4+2*y2^2*x3^3*r1^4*x2+2*y2^2*x3^3*r1^2*x2^3-4*y2^2*x3^4*r1^2*x2^2+2*y2^2*x3^4*r1^2*r2^2+2*y2^2*x3^5*r1^2*x2+2*y2^2*x3^4*x2^2*r2^2-2*y2^2*x3^3*x2^3*r3^2-2*y2^2*x3^5*r2^2*x2-3*y2^2*x3^4*x2^2*y3^2+2*y2^2*x3^4*x2^2*r3^2+2*x2^4*y3^3*y2*r1^2-2*x2^4*y3^3*y2*x3^2-2*x2^4*y3^3*y2*r3^2+2*x2^3*y3^4*r1^2*x3-2*x2^3*y3^4*y2^2*x3-2*x2^3*y3^4*r2^2*x3-2*x2^2*y3^3*y2^3*r3^2+2*x2^2*y3^4*y2^2*r2^2+2*x2^2*y3^5*r1^2*y2+2*x2^2*y3^3*r1^4*y2-4*x2^2*y3^4*r1^2*y2^2+2*x2^2*y3^4*r1^2*r2^2+2*x2^2*y3^3*y2^3*r1^2+2*x2^2*y3^4*y2^2*r3^2-2*x2^2*y3^5*y2*r2^2-x2^2*y3^2*y2^2*r1^4-x2^2*y3^2*y2^2*r3^4-x2^2*y3^2*r1^4*x3^2-x2^2*y3^2*r2^4*x3^2+2*x2^3*y3^2*r1^4*x3+2*x2^5*y3^2*r1^2*x3+2*x2^2*y3^2*r1^2*x3^2*r2^2-8*x2^3*y3^3*r1^2*y2*x3+2*y3^2*z2^2*r1^2*x2^2*r3^2+2*y3^2*z2^2*x2^2*z3^2*r3^2+4*y2^4*x3*x2*y3^2*r3^2+4*y2^3*x3*x2*y3^3*z2^2-4*y2^3*x3*x2*y3^3*r2^2-4*y2^2*x3*x2*y3^4*r1^2-4*y2^2*x3*x2*y3^2*r1^4+8*y2^3*x3*x2*y3^3*r1^2+4*y2*x3*x2^3*y3^3*z2^2-4*y2*x3*x2^3*y3^3*r2^2-4*y2^4*x3*x2*y3^2*r1^2+4*y2^3*x3^3*x2*y3*z3^2-4*y2^3*x3^3*x2*y3*r3^2+4*y2^3*x3*x2*y3^3*z3^2-4*y2^3*x3*x2*y3^3*r3^2+4*y2^2*x3*x2*y3^4*r2^2+2*y2*x3*x2*y3^3*r1^4+2*y2^3*x3*x2*y3*r1^4+2*y2^3*x3*x2*y3*z3^4+2*y2^3*x3*x2*y3*r3^4+2*y2*x3*x2*y3^3*z2^4+2*y2*x3*x2*y3^3*r2^4+2*y2*x3*x2^3*y3*r1^4+2*y2*x3^3*x2*y3*r1^4+2*y2*x3^3*x2*y3*z2^4+2*y2*x3^3*x2*y3*r2^4+2*y2*x3*x2^3*y3*z3^4+2*y2*x3*x2^3*y3*r3^4-4*y2*x3^2*x2^2*y3*r1^4-4*y2*x3^2*x2^4*y3*r1^2+8*y2*x3^3*x2^3*y3*r1^2-4*y2*x3^4*x2^2*y3*r1^2+4*y2*x3^3*x2^3*y3*z2^2-4*y2*x3^3*x2^3*y3*r2^2+4*y2*x3^2*x2^4*y3*r3^2+4*y2^3*x3^3*x2*y3*z2^2-4*y2^3*x3^3*x2*y3*r2^2+4*y2*x3^4*x2^2*y3*r2^2+4*y2*x3^3*x2^3*y3*z3^2-4*y2*x3^3*x2^3*y3*r3^2+4*y2*x3*x2^3*y3^3*z3^2-4*y2*x3*x2^3*y3^3*r3^2+16*y2^2*x3^2*x2^2*y3^2*r1^2);
b=(-z2^3*y3^2-x2^2*y3^2*z2-y2^2*z3*x3^2-y2^2*z3*y3^2+y2^3*z3*y3+y2*y3^3*z2-y2^2*y3^2*z2-z2*x3^2*x2^2-z2*x3^2*y2^2+z2*x3^3*x2+x2^3*z3*x3-x2^2*z3*x3^2-x2^2*z3*y3^2+y2*z3*z2^2*y3+y2*x3^2*z2*y3+y2*z3^2*z2*y3+z2*x3*x2*y3^2+z2*x3*x2*z3^2+x2*z3*y2^2*x3+x2*z3*z2^2*x3+x2^2*y3*y2*z3-y2^2*z3^3-z2^3*x3^2-x2^2*z3^3-r1^2*y3^2*z2-y2^2*z3*r1^2+r1^2*y3*y2*z3+y2*r1^2*z2*y3+z2*y3^2*r2^2-z2*x3^2*r1^2+z2*x3^2*r2^2-x2^2*z3*r1^2+x2^2*z3*r3^2+y2^2*z3*r3^2-y2*z3*r2^2*y3-y2*r3^2*z2*y3+z2*x3*r1^2*x2-z2*x3*x2*r3^2+x2*z3*r1^2*x3-x2*z3*r2^2*x3);
c=(-2*y2*z3*z2*y3-2*z2*x3*x2*z3+z3^2*y2^2+x3^2*z2^2+z3^2*x2^2+y2^2*x3^2+x2^2*y3^2+y3^2*z2^2-2*y2*x3*x2*y3);
if(a<0||c==0), result=[nan;nan;nan]; return; end % coz c is the denominator and a is under a root
%     error('Error in interx.m at z'); end
za=-1/2*(b-a^(1/2))/c;
zb=-1/2*(b+a^(1/2))/c;
if(za>zb)
    if(pos), z=za; else z=zb; end
else
    if(pos), z=zb; else z=za; end
end
a=(2*z*z2*x3-2*x2*z*z3+r1^2*x2-r1^2*x3-x2^2*x3-y2^2*x3-z2^2*x3+r2^2*x3+x2*x3^2+x2*y3^2+x2*z3^2-x2*r3^2);
b=(-2*y2*x3+2*x2*y3);
if(b==0), result=[nan;nan;nan]; return; end % coz b is the denominator in the expression
%     error('Error in interx.m at y'); end
y=a/b;
if(x2==0), result=[nan;nan;nan]; return; end
%     error('Error in interx.m at x'); end
x = 1/2*(r1^2+x2^2-2*y*y2+y2^2-2*z*z2+z2^2-r2^2)/x2;
%% convert result back to global
result=[x1;y1;z1;1]+[x;y;z;0];
% disp([x1 y1 z1 x2+x1 y2+y1 z2+z1 x3+x1 y3+y1 z3+z1]);
%disp('Solution'); disp(result');
%disp('Constraints'); % check distances to three centers
%disp((result(1)-X1(1))^2+(result(2)-X1(2))^2+(result(3)-X1(3))^2-r1^2);
%disp((result(1)-X2(1))^2+(result(2)-X2(2))^2+(result(3)-X2(3))^2-r2^2);
%disp((result(1)-X3(1))^2+(result(2)-X3(2))^2+(result(3)-X3(3))^2-r3^2);