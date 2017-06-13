load output_flasso.mat
E_min   = min(E_PD3O_min); % computed from 10,000 iterations using PD3O
addpath('fcns','data','output') 
plot_res = 'energy';
plot_res = 'LS';
plot_res = 'PD';


switch plot_res 
    case 'energy'
        cs={'r-','b-','k-','m-','r--','b--','k--','m--','r-.','b-.','k-.','m-.'};
        h1  = figure;
        legend_text = {};

        j = 1;
        figure(h1)
        semilogy(E_PD3O(:)./E_min-1,cs{j});
        ylim([1e-14,1e3])
        xlim([1,4e3])
        legend_text(end+1) = {'PD3O-$\gamma_1$'};

        j = 2;
        figure(h1)
        hold on
        semilogy(E_PDFP(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_1$'};

        j = 3;
        figure(h1)
        hold on
        semilogy(E_CV(:)./E_min-1,cs{j});
        legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};

        j = 4;
        figure(h1)
        hold on
        semilogy(E_AFBA(:)./E_min-1,cs{j});
        legend_text(end+1) = {'AFBA-$\gamma_1$'};

        j = 5;
        figure(h1)
        hold on
        semilogy(E_PD3O2(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PD3O-$\gamma_2$'};

        j = 6;
        figure(h1)
        hold on
        semilogy(E_PDFP2(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_2$'};

%         j = 7;
%         figure(h1)
%         hold on
%         semilogy(E_CV2(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'Condat-Vu-$\gamma_2$'};
% 
%         j = 8;
%         figure(h1)
%         hold on
%         semilogy(E_AFBA2(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'AFBA-$\gamma_2$'};


        j = 9;
        figure(h1)
        hold on
        semilogy(E_PD3O3(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PD3O-$\gamma_3$'};

        j = 10;
        figure(h1)
        hold on
        semilogy(E_PDFP3(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_3$'};

%         j = 11;
%         figure(h1)
%         hold on
%         semilogy(E_CV3(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'Condat-Vu-$\gamma_3$'};
% 
%         j = 12;
%         figure(h1)
%         hold on
%         semilogy(E_AFBA3(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'AFBA-$\gamma_3$'};

        h_legend = legend(legend_text,'Interpreter','latex');
        set(h_legend,'FontSize',10);
        xlabel('iteration','FontSize',20)
        ylabel('$\frac{f-f^*}{f^*}$','Interpreter','LaTex','FontSize',20);

        myprint('output/flasso_energy_1',h1)   % print the file in .PDf and .eps
        
        h2  = figure;
        legend_text = {};
        j = 1;
        figure(h2)
        semilogy(E_PD3O4(:)./E_min-1,cs{j});
        ylim([1e-14,1e3])
        xlim([1,4e3])
        legend_text(end+1) = {'PD3O-$\lambda_1$'};

        j = 2;
        figure(h2)
        hold on
        semilogy(E_PDFP4(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_1$'};

        j = 3;
        figure(h2)
        hold on
        semilogy(E_CV4(:)./E_min-1,cs{j});
        legend_text(end+1) = {'Condat-Vu-$\lambda_1$'};

        j = 4;
        figure(h2)
        hold on
        semilogy(E_AFBA4(:)./E_min-1,cs{j});
        legend_text(end+1) = {'AFBA-$\lambda_1$'};
        
        j = 5;
        figure(h2)
        hold on
        semilogy(E_PD3O5(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PD3O-$\lambda_2$'};

        j = 6;
        figure(h2)
        hold on
        semilogy(E_PDFP5(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_2$'};

%         j = 7;
%         figure(h2)
%         hold on
%         semilogy(E_CV5(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'Condat-Vu-$\lambda_2$'};
% 
%         j = 8;
%         figure(h2)
%         hold on
%         semilogy(E_AFBA5(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'AFBA-$\lambda_2$'};
        
        j = 9;
        figure(h2)
        hold on
        semilogy(E_PD3O6(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PD3O-$\lambda_3$'};

        j = 10;
        figure(h2)
        hold on
        semilogy(E_PDFP6(:)./E_min-1,cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_3$'};

%         j = 11;
%         figure(h2)
%         hold on
%         semilogy(E_CV6(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'Condat-Vu-$\lambda_3$'};
% 
%         j = 12;
%         figure(h2)
%         hold on
%         semilogy(E_AFBA6(:)./E_min-1,cs{j});
%         legend_text(end+1) = {'AFBA-$\lambda_3$'};
        

        h_legend = legend(legend_text,'Interpreter','latex');
        set(h_legend,'FontSize',10);
        xlabel('iteration','FontSize',20)
        ylabel('$\frac{f-f^*}{f^*}$','Interpreter','LaTex','FontSize',20);

        myprint('output/flasso_energy_2',h1)   % print the file in .PDf and .eps        
        
    case 'LS'
        cs={'r-','b-','k-','m-','r--','b--','k--','m--','r-.','b-.','k-.','m-.'};
        h1  = figure;
        legend_text = {};

        j = 1;
        figure(h1)
        semilogy(out_PD3O.LS(:),cs{j});
        ylim([1e-12,1e3])
        xlim([1,4e3])
        legend_text(end+1) = {'PD3O-$\gamma_1$'};

        j = 2;
        figure(h1)
        hold on
        semilogy(out_PDFP.LS(:),cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_1$'};

        j = 3;
        figure(h1)
        hold on
        semilogy(out_CV.LS(:),cs{j});
        legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};

        j = 4;
        figure(h1)
        hold on
        semilogy(out_AFBA.LS(:),cs{j});
        legend_text(end+1) = {'AFBA-$\gamma_1$'};
        
        j = 5;
        figure(h1)
        hold on
        semilogy(out_PD3O2.LS(:),cs{j});
        legend_text(end+1) = {'PD3O-$\gamma_2$'};

        j = 6;
        figure(h1)
        hold on
        semilogy(out_PDFP2.LS(:),cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_2$'};

        j = 9;
        figure(h1)
        hold on
        semilogy(out_PD3O3.LS(:),cs{j});
        legend_text(end+1) = {'PD3O-$\gamma_3$'};

        j = 10;
        figure(h1)
        hold on
        semilogy(out_PDFP3.LS(:),cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_3$'};

        h_legend = legend(legend_text,'Interpreter','latex');
        set(h_legend,'FontSize',10);
        xlabel('iteration','FontSize',20)
        ylabel('$\|x-x^*\|^2$','Interpreter','LaTex','FontSize',20);

        myprint('output/flasso_LS_1',h1)   % print the file in .PDf and .eps        
        
        h2  = figure;
        legend_text = {};

        j = 1;
        figure(h2)
        semilogy(out_PD3O4.LS(:),cs{j});
        ylim([1e-12,1e3])
        xlim([1,4e3])
        legend_text(end+1) = {'PD3O-$\lambda_1$'};

        j = 2;
        figure(h2)
        hold on
        semilogy(out_PDFP4.LS(:),cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_1$'};

        j = 3;
        figure(h2)
        hold on
        semilogy(out_CV4.LS(:),cs{j});
        legend_text(end+1) = {'Condat-Vu-$\lambda_1$'};

        j = 4;
        figure(h2)
        hold on
        semilogy(out_AFBA4.LS(:),cs{j});
        legend_text(end+1) = {'AFBA-$\lambda_1$'};
        
        j = 5;
        figure(h2)
        hold on
        semilogy(out_PD3O5.LS(:),cs{j});
        legend_text(end+1) = {'PD3O-$\lambda_2$'};

        j = 6;
        figure(h2)
        hold on
        semilogy(out_PDFP5.LS(:),cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_2$'};

        j = 9;
        figure(h2)
        hold on
        semilogy(out_PD3O6.LS(:),cs{j});
        legend_text(end+1) = {'PD3O-$\lambda_3$'};

        j = 10;
        figure(h2)
        hold on
        semilogy(out_PDFP6.LS(:),cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_3$'};

        h_legend = legend(legend_text,'Interpreter','latex');
        set(h_legend,'FontSize',10);
        xlabel('iteration','FontSize',20)
        ylabel('$\|x-x^*\|^2$','Interpreter','LaTex','FontSize',20);

        myprint('output/flasso_LS_2',h1)   % print the file in .PDf and .eps               
    case 'PD'
        cs={'r-','b-','k-','m-','r--','b--','k--','m--','r-.','b-.','k-.','m-.'};
        h1  = figure;
        legend_text = {};

        j = 1;
        figure(h1)
        semilogy(out_PD3O.PD(:),cs{j});
        xlim([1,4e3])
        ylim([1e-9,1e6])
        legend_text(end+1) = {'PD3O-$\gamma_1$'};

        j = 2;
        figure(h1)
        hold on
        semilogy(out_PDFP.PD(:),cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_1$'};

        j = 3;
        figure(h1)
        hold on
        semilogy(out_CV.PD(:),cs{j});
        legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};

        j = 4;
        figure(h1)
        hold on
        semilogy(out_AFBA.PD(:),cs{j});
        legend_text(end+1) = {'AFBA-$\gamma_1$'};

        j = 5;
        figure(h1)
        hold on
        semilogy(out_PD3O2.PD(:),cs{j});
        legend_text(end+1) = {'PD3O-$\gamma_2$'};

        j = 6;
        figure(h1)
        hold on
        semilogy(out_PDFP2.PD(:),cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_2$'};

        j = 9;
        figure(h1)
        hold on
        semilogy(out_PD3O3.PD(:),cs{j});
        legend_text(end+1) = {'PD3O-$\gamma_3$'};

        j = 10;
        figure(h1)
        hold on
        semilogy(out_PDFP3.PD(:),cs{j});
        legend_text(end+1) = {'PDFP-$\gamma_3$'};

        h_legend = legend(legend_text,'Interpreter','latex');
        set(h_legend,'FontSize',10);
        xlabel('iteration','FontSize',20)
        ylabel('primal-dual gap','Interpreter','LaTex','FontSize',20);

        myprint('output/flasso_PD_1',h1)   % print the file in .PDf and .eps        
        h2  = figure;
        legend_text = {};

        j = 1;
        figure(h2)
        semilogy(out_PD3O4.PD(:),cs{j});
        xlim([1,4e3])
        ylim([1e-9,1e6])
        legend_text(end+1) = {'PD3O-$\lambda_1$'};

        j = 2;
        figure(h2)
        hold on
        semilogy(out_PDFP4.PD(:),cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_1$'};

        j = 3;
        figure(h2)
        hold on
        semilogy(out_CV4.PD(:),cs{j});
        legend_text(end+1) = {'Condat-Vu-$\lambda_1$'};

        j = 4;
        figure(h2)
        hold on
        semilogy(out_AFBA4.PD(:),cs{j});
        legend_text(end+1) = {'AFBA-$\lambda_1$'};

        j = 5;
        figure(h2)
        hold on
        semilogy(out_PD3O5.PD(:),cs{j});
        legend_text(end+1) = {'PD3O-$\lambda_2$'};

        j = 6;
        figure(h2)
        hold on
        semilogy(out_PDFP5.PD(:),cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_2$'};

        j = 9;
        figure(h2)
        hold on
        semilogy(out_PD3O6.PD(:),cs{j});
        legend_text(end+1) = {'PD3O-$\lambda_3$'};

        j = 10;
        figure(h2)
        hold on
        semilogy(out_PDFP6.PD(:),cs{j});
        legend_text(end+1) = {'PDFP-$\lambda_3$'};

        h_legend = legend(legend_text,'Interpreter','latex');
        set(h_legend,'FontSize',10);
        xlabel('iteration','FontSize',20)
        ylabel('primal-dual gap','Interpreter','LaTex','FontSize',20);

        myprint('output/flasso_PD_2',h1)   % print the file in .PDf and .eps             
end