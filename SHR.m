
% This program solves MHD stagnation point flow of nanofluid
% f'''+ff''-f'^2 + M*(\lambda-f') + lmbda^2=0
% theta''+ pr(Nbphi'theta'+Nt theta'^2)=0
% phi''+Lefphi'+Nt/Nb theta''=0
% with boundary conditions
%f(0)=0,f'(0)=1,theta(0)=1,phi(0)=0 at eta=0
%[y0(1);y0(2)-1;y0(4)-1;y0(6)-1;yinf(2)-A;yinf(4);yinf(6)];
%f'(\infinity)=lambda, \theta(\infinity)=0, \phi(\infinity)=0 as eta tends to
% infinity
 function [solu] = SHR

global A Pr Nb Nt Le M S L 

infinity = 5;

% Initial values
solu1 = bvpinit(linspace(1e-7,infinity,100),zeros(7,1));
%solu1 = bvpinit(linspace(0,10,5),zeros(7,1));
% solinit=bvpinit(linspace(0,infinity,8),[0 1 1]);
% bvp4c takes two functions defined below and give solution in structure
% form
 options = bvpset('RelTol',1e-9,'AbsTol',1e-9);
 solu = bvp5c(@bvpex1, @bcex1,solu1,options);

 x = solu.x;
 
 y = solu.y;

plot(x,y(2,:),'b')
%xlabel('\eta')
%ylabel('df/dx')

%hold on;

% find all y(0), y'(0), y''(0)
%deval(sol,0:6)
%y(3,1)

% Here I defined first order ODEs 

function ysolu = bvpex1(x,y)
  ysolu = [y(2);
          y(3);
           - y(1) * y(3) + y(2)*y(2)  + M * y(2)-M*A-A*A;
           y(5);
           - Pr * (y(1)*y(5)  + Nb * y(5) * y(7)+ Nt * y(5)*y(5));
           y(7);
           - Le * y(1) * y(7) + Nt* Pr * (y(1) * y(5) + Nb * y(5) * y(7)+ Nt * y(5)*y(5))/Nb];
%     end
end

% residual of boundary conditions
function res = bcex1(y0, yinf)

res = [y0(1)-S;
       y0(2)-L;
       y0(4)-1;
       y0(6)-1;
       yinf(2)-A;
       yinf(4);
       yinf(6)];
end
%fprintf('i=% d',i);
 % plot(M,y(2,:))
%  %plot(x, y(2, :),'g')
%  xlabel('m')
% display x
end