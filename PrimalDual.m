classdef PrimalDual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Primal-Dual algorithms for minimizing                            
%    F(x) + G(x) + H(Ax)                                              
%    F is differentiable, and both G and H are proximable             
%
%
%    Three methods are implemented: Condat-Vu, PDFP, and PD3O        
%    Reference: 
%       M. Yan, "A Primal-Dual Three Operator Splitting Scheme", arxiv
%       1611.????, Nov., 2016.
%
%    Contact:
%       Ming Yan yanm @ math.msu.edu
%       Downloadable from https://github.com/mingyan08/PD3O                                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   properties
      lambda;       % parameter lambda
      gamma;        % parameter gamma
      myF       =  @(x)   0;            % smooth function    F: x --> F(x)
      myG       =  @(x)   0;            % proximable fcn     G: x --> G(x)
      myH       =  @(y)   0;            % proximable fcn     H: y --> H(y)
      myA       =  @(x)   0;            % linear operator    A: x --> A(x)
      myGradF   =  @(x)   0;            % gradient of F      gradF: x    --> grad(F)(x)
      myProxG   =  @(x,t) x;            % prox of G          proxG: x,t  --> prox(t.G)(x)
      myProxH   =  @(y,t) y;            % prox of H          proxH: y,t  --> prox(t.H)(y)
      myAdjA    =  @(y)   0;            % adjoint of A       adjA:  y    --> A'*(y)
   end
   methods
       function s = E(this, x)  
           % return the function value for given x: F(x) + G(x) + H(Ax) 
           s = this.myF(x) + this.myG(x) + this.myH(this.myA(x));
       end
       
       function [x, s, E] = minimize(this, x, s, iter, method)       
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % s^+ = (I - Prox_{gamma/lambda H})(s + A * x_bar)
           % x^+ = Prox_{gamma G}(x - gamma GradF - lambda A' * s^+)
           %     = Prox_{gramm G}(xh - lambda A' * s^+)
           % x_bar^+ = 2 x^+ - x + gamma GradF - gamma GradF'
           %         = 2 x^+ - xh - gamma GradF'
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           x_bar = x;
           E     = zeros(iter, 1);
           gradF = this.myGradF(x);
           mu    = this.gamma/this.lambda;
           for i = 1:iter
               x0       = x;
               xh       = x - this.gamma * gradF;  % x - \gamma \nabla F(x)
               sh       = s + this.myA(x_bar);     % s + A * x_bar
               s        = sh - this.myProxH(sh, mu);
               x        = this.myProxG(xh - this.lambda * this.myAdjA(s), this.gamma);
               gradF    = this.myGradF(x);
               switch method
                   case 'PD3O'  % Primal-Dual Three Operator
                       x_bar    = 2 * x - xh - this.gamma * gradF;
                   case {'CV','Condat-Vu'}  % Condat-Vu 
                       x_bar    = 2 * x - x0;
                   case 'PDFP'  % Primal Dual Fixed Point 
                       x_bar    = this.myProxG(x - this.gamma * gradF - this.lambda * this.myAdjA(s), this.gamma);
                   otherwise
                       warning('Unexpected Method, Please choose from PDFP, PD3O, and CP')
               end
               E(i)     = this.E(x);
           end
       end
   end
end