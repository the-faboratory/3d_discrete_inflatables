function [T,N,B,k,t] = frenet(x,y,z)
% FRENET - Frenet-Serret Space Curve Invarients
%
%   [T,N,B,k,t] = frenet(x,y);
%   [T,N,B,k,t] = frenet(x,y,z);
%
%   Returns the 3 vector and 2 scaler invarients of a space curve defined
%   by vectors x,y and z.  If z is omitted then the curve is only a 2D,
%   but the equations are still valid.
%
%    _    r'
%    T = ----  (Tangent)
%        |r'|
%
%    _    T'
%    N = ----  (Normal)
%        |T'|
%    _   _   _
%    B = T x N (Binormal)
%
%    k = |T'|  (Curvature)
%
%    t = vdot(cross(dr, ddr), dddr) ./ mag(cross(dr, ddr),1).^2; (torsion)
%
%
%    Example:
%    theta = 2*pi*linspace(0,2,100);
%    x = cos(theta);
%    y = sin(theta);
%    z = theta/(2*pi);
%    [T,N,B,k,t] = frenet(x,y,z);
%    line(x,y,z); hold on
%    quiver3(x,y,z,T(:,1)',T(:,2)',T(:,3)','color','r')
%    quiver3(x,y,z,N(:,1)',N(:,2)',N(:,3)','color','g')
%    quiver3(x,y,z,B(:,1)',B(:,2)',B(:,3)','color','b')
%    legend('Curve','Tangent','Normal','Binormal')
%
%
% See also: GRADIENT

if nargin == 2
    z = zeros(size(x));
end

% CONVERT TO COLUMN VECTOR
x = x(:);
y = y(:);
z = z(:);

% SPEED OF CURVE
dx = gradient(x);
dy = gradient(y);
dz = gradient(z);
dr = [dx dy dz];

% Accelleration OF CURVE
ddx = gradient(dx);
ddy = gradient(dy);
ddz = gradient(dz);
ddr = [ddx ddy ddz];

dddx = gradient(ddx);
dddy = gradient(ddy);
dddz = gradient(ddz);
dddr = [dddx dddy dddz];

% TANGENT
T = dr./mag(dr,3);

%ROB ADD:process T vectors: for noise as well. In particular we look for cases
%where components of T, N, or B are quite close to zero (i.e. 0.000001 or
%-0.000001) based on rounding, and just map them to be all the same sign.
%Normally TNB is ill defined for planar curves, but here we abide by the
%special case as articulated in: https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas
for i = 1:length(T)
    for j = 1:3
        if abs(T(i,j)) < 1e-5
            T(i,j) = abs(T(i,j));
        end
    end
end
%%%%%END ROB ADD

% DERIVIATIVE OF TANGENT
dTx =  gradient(T(:,1));
dTy =  gradient(T(:,2));
dTz =  gradient(T(:,3));

dT = [dTx dTy dTz];


% NORMAL
N = dT./mag(dT,3);
% BINORMAL
B = cross(T,N);
% CURVATURE
% k = mag(dT,1);
k = mag(cross(dr,ddr),1)./((mag(dr,1)).^3);  % ROB: validated equation here https://math.libretexts.org/Bookshelves/Calculus/Supplemental_Modules_(Calculus)/Vector_Calculus/2%3A_Vector-Valued_Functions_and_Motion_in_Space/2.3%3A_Curvature_and_Normal_Vectors_of_a_Curve
% TORSION
%size(mag(cross(dr, ddr),1).^2)
%size(vdot(cross(dr, ddr), dddr))
t = vdot(cross(dr, ddr), dddr) ./ mag(cross(dr, ddr),1).^2;  % ROB: and here https://www.math.upenn.edu/~wziller/math114f13/ch13-5+6.pdf
%t = vdot([gradient(N(:,1)) gradient(N(:,2)) gradient(N(:,3)) ],B); % ROB: can also express this way. 

function N = vdot(A, B)
%row-wise dot-product of A and B
N=zeros(size(A,1),1);
for i=1:size(A,1)
    N(i) = dot(A(i,:), B(i,:));
end


function N = mag(T,n)
% MAGNATUDE OF A VECTOR (Nx3)
%  M = mag(U)
N = sum(abs(T).^2,2).^(1/2);
d = find(N==0);
N(d) = eps*ones(size(d));
N = N(:,ones(n,1));