% Post-processing of this run folder: energy with contractility recruitment
% marked, contractility and purse-string tension vs time, per-step CSV.
folder=fileparts(mfilename('fullpath'));
s=load(fullfile(folder,'workspace.mat'));
n=s.param.Nsteps;steps=(1:n)';t=s.timemat(:);c=s.contractilityFactorLog(:);e=s.energymat(:);
dc=diff([1;c]);increased=dc>1e-12;
assert(numel(c)==n && all(isfinite(c)) && all(increased));
% Reconstruct the purse-string tension from the exact original recurrence.
lambda=zeros(n,1);previous=0;
for k=1:n
 recruited=1-exp(-k*s.param.deltat/s.param.tau_Kw);
 previous=previous+s.param.deltat*s.param.k_recruit*(s.param.lambda_purse_string*recruited-previous);
 lambda(k)=previous;
end
assert(abs(lambda(end)-s.celldata.lambda_current)<1e-10);
assert(max(abs(c-(1+(s.param.contractility_wound-1)*lambda/s.param.lambda_purse_string)))<1e-10);
sample=(100:100:n)';
set(groot,'defaultFigureColor','white','defaultAxesColor','white',...
 'defaultAxesXColor','black','defaultAxesYColor','black','defaultTextColor','black');
blue=[0.05 0.32 0.68];red=[0.78 0.16 0.12];
f=figure('Visible','off','Position',[60 60 1200 800]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(steps,e,'Color',blue,'LineWidth',1.6);hold on;
scatter(sample,e(sample),32,red,'filled');
for k=sample';xline(k,':','Color',[0.8 0.8 0.8],'HandleVisibility','off');end
grid on;xlim([0 n]);ylabel('Mechanical potential energy');
title('Energy vs timestep: contractility recruitment marked','Color','black');
legend({'Saved energy','Steps sampled every 100 (not discrete jumps)'},'Location','northeast','TextColor','black','Color','white');
text(0.02,0.07,'Contractility increases at every step, 1-1000. Red points mark every 100th step.',...
 'Units','normalized','FontSize',10,'Color',red);
nexttile;
yyaxis left;plot([0;steps],[1;c],'Color',red,'LineWidth',1.8);hold on;
scatter(sample,c(sample),28,red,'filled');ylabel('Margin perimeter-contractility factor (x bulk)');
ax=gca;ax.YAxis(1).Color=red;
yyaxis right;plot(steps,dc,'Color',[0.3 0.3 0.3],'LineWidth',1.2);ylabel('Increase per timestep');ax.YAxis(2).Color=[0.3 0.3 0.3];
xlim([0 n]);grid on;xlabel('Timestep (one step = 0.01 simulation-time units)');
title('Continuous recruitment: increasing level with a progressively smaller increment','Color','black');
f.Color='white';
exportgraphics(f,fullfile(folder,'energy_vs_time_contractility_marked.png'),'Resolution',200,'BackgroundColor','white');
savefig(f,fullfile(folder,'energy_vs_time_contractility_marked.fig'));close(f);

f=figure('Visible','off','Position',[60 60 1200 800]);tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;plot([0;t],[1;c],'Color',red,'LineWidth',2);hold on;
scatter(t(sample),c(sample),28,red,'filled');yline(s.param.contractility_wound,'--','Target = 3x','Color',[0.35 0.35 0.35]);
grid on;xlim([0 n*s.param.deltat]);ylim([0.9 3.15]);ylabel('Margin perimeter-contractility factor (x bulk)');
xlabel('Simulation time');title('Wound-margin contractility vs time','Color','black');
text(0.45,1.2,sprintf('Initial: 1.000x   |   final: %.6fx\nIncreases at all 1000 timesteps',c(end)),'Color','black');
nexttile;plot([0;t],[0;lambda],'Color',blue,'LineWidth',2);hold on;
yline(s.param.lambda_purse_string,'--','Target = 4','Color',[0.35 0.35 0.35]);
grid on;xlim([0 n*s.param.deltat]);ylim([0 4.2]);xlabel('Simulation time');ylabel('Purse-string line tension');
title('Purse-string tension (separate from the cell contractility factor)','Color','black');
f.Color='white';
exportgraphics(f,fullfile(folder,'contractility_vs_time.png'),'Resolution',200,'BackgroundColor','white');
savefig(f,fullfile(folder,'contractility_vs_time.fig'));close(f);

T=table(steps,t,e,c,dc,increased,c/s.param.rstiff,lambda,...
 'VariableNames',{'timestep','simulation_time','energy','margin_contractility_factor',...
 'contractility_increase_this_step','contractility_increased','perimeter_energy_coefficient','purse_string_tension_reconstructed'});
writetable(T,fullfile(folder,'contractility_vs_time_all_steps.csv'));
fid=fopen(fullfile(folder,'contractility_plot_notes.txt'),'w');
fprintf(fid,'Computed from workspace.mat of this run.\n');
fprintf(fid,'Margin perimeter-contractility factor (contractilityFactorLog) rises at every timestep 1-1000 (t = 0.01 to 10), from 1 to %.12f.\n',c(end));
fprintf(fid,'Energy plot: red points mark every 100th step; per-step values are in the CSV.\n');
fprintf(fid,'Purse-string tension from the recruitment recurrence (matches the saved final value): final tension %.12f.\n',lambda(end));
fprintf(fid,'Perimeter energy coefficient = contractility factor / rstiff: bulk 2, final margin %.12f (model units).\n',c(end)/s.param.rstiff);
fclose(fid);
fprintf('Increases: %d/%d; first step=%d; last step=%d; final contractility=%.12f; final lambda=%.12f\n',nnz(increased),n,find(increased,1),find(increased,1,'last'),c(end),lambda(end));
fprintf('Saved annotated energy, contractility plots, editable FIGs, per-step CSV and notes in %s\n',folder);
