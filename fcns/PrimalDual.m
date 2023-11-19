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
      input;
      myF       =  @(x)   0;            % smooth function    F: x --> F(x)
      myF2      =  @(x)   0;
      myG       =  @(x)   0;            % proximable fcn     G: x --> G(x)
      myH       =  @(y)   0;            % proximable fcn     H: y --> H(y)
      myHA      =  @(y)   0;            % conjugate fcn of H
      myA       =  @(x)   0;            % linear operator    A: x --> A(x)
      myGradF   =  @(x)   0;            % gradient of F      gradF: x    --> grad(F)(x)
      myProxG   =  @(x,t) x;            % prox of G          proxG: x,t  --> prox(t.G)(x)
      myProxH   =  @(y,t) y;            % prox of H          proxH: y,t  --> prox(t.H)(y)
      myAdjA    =  @(y)   0;            % adjoint of A       adjA:  y    --> A'*(y)
   end
   methods
       function p = E(this, x)  
           % return the function value for given x: F(x) + G(x) + H(Ax) 
           p = this.myF(x) + this.myG(x) + this.myH(this.myA(x));
       end
       
       function pd = gap2(this, x, s, x_min, s_min)
           temp = this.lambda/this.gamma;
           % return the primal dual gap: F(x) + G(x) + <s, Ax> - H*(s)
           pd1 = this.myF2(x) + this.myG(x) + trace(s_min'*this.myA(x)) - this.myHA(s_min); 
           pd2 = this.myF2(x_min) + this.myG(x_min) + temp*trace(s'*this.myA(x_min)) - this.myHA(temp*s);
           pd  = pd1 - pd2;           
       end
       
       function pd = gap(this, x, s, x_min, s_min)
           temp = this.lambda/this.gamma;
           % return the primal dual gap: F(x) + G(x) + <s, Ax> - H*(s)
           pd1 = this.myF(x) + this.myG(x) + trace(s_min'*this.myA(x)) - this.myHA(s_min); 
           pd2 = this.myF(x_min) + this.myG(x_min) + temp*trace(s'*this.myA(x_min)) - this.myHA(temp*s);
           pd  = pd1 - pd2;           
       end
       
       function [x, s, E, out] = minimize(this, method, type)       
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % s^+ = (I - Prox_{gamma/lambda H})(s + A * x_bar)
           % x^+ = Prox_{gamma G}(x - gamma GradF - lambda A' * s^+)
           %     = Prox_{gramm G}(xh - lambda A' * s^+)
           % x_bar^+ = 2 x^+ - x + gamma GradF - gamma GradF'
           %         = 2 x^+ - xh - gamma GradF'
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if ~isfield(this.input, 'x')
                error('Error. Put the initial primal x in this.input.x')               
           end
           x     = this.input.x;
           if ~isfield(this.input, 's')
                error('Error. Put the initial dual s in this.input.s')               
           end
           s     = this.input.s;
           if isfield(this.input, 'iter')
               iter  = this.input.iter;
           else
               iter = 1000;
           end
           if isfield(this.input, 'x_min')
               x_min = this.input.x_min;
           end
           if isfield(this.input, 's_min')
               s_min = this.input.s_min;
           end
           
           x_bar = x;
           E     = zeros(iter, 1);
           gradF = this.myGradF(x);
           mu    = this.gamma/this.lambda;
           out   = struct();

           compare = 1;
           if (~exist('x_min', 'var') || isempty(x_min))
               compare = 0;
           end;
           if exist('x_min', 'var') || ~isempty(type)
               if type == 1;
                   out.LS = zeros(iter,1);
               end
               pdgap = 0;
               if exist('s_min', 'var')
                   pdgap = 1;
                   out.PD = zeros(iter,1);
                   out.PD2 = zeros(iter,1);
               end
           end
           
           
           for i = 1:iter
               x0       = x;
               s0       = s;
               xh       = x - this.gamma * gradF;  % x - \gamma \nabla F(x)
               sh       = s + this.myA(x_bar);     % s + A * x_bar
               s        = sh - this.myProxH(sh, mu);
               if strcmp(method,'AFBA')
                   x        = x_bar - this.lambda * this.myAdjA(s - s0);
               else
                   x        = this.myProxG(xh - this.lambda * this.myAdjA(s), this.gamma);
               end
               gradF    = this.myGradF(x);
               switch method
                   case 'PD3O'  % Primal-Dual Three Operator
                       x_bar    = 2 * x - xh - this.gamma * gradF;
                   case {'CV','Condat-Vu'}  % Condat-Vu 
                       x_bar    = 2 * x - x0;
                   case {'PDFP', 'AFBA'}  % Primal Dual Fixed Point and AFBA
                       x_bar    = this.myProxG(x - this.gamma * gradF - this.lambda * this.myAdjA(s), this.gamma);
                   otherwise
                       warning('Unexpected Method, Please choose from PDFP, PD3O, and CP ')
               end
               E(i)     = this.E(x);
               if compare 
                   switch type 
                       case 1
                           out.LS(i) = norm(x - x_min, 'fro');
                           if pdgap
                               out.PD(i) = this.gap(x0, s, x_min, s_min);
                               out.PD2(i) = this.gap2(x0, s, x_min, s_min);                               
                           end
                   end
               end
           end
           s = (this.lambda/this.gamma)*s;
       end
   end
end
