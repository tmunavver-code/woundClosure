folder='/Users/tmunavver/Desktop/woundclosure/2D_VX_model/results/force_closure_gradual_intercalation_crawling_active_1000_steps_2026-10-24_002857';
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
text(0.02,0.07,'Contractility increases at EVERY step, 1–1000. Red points label a readable subset.',...
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
fprintf(fid,'Derived from this folder''s original workspace.mat without running the simulation or altering original outputs.\n');
fprintf(fid,'Cell perimeter-contractility factor is directly saved in contractilityFactorLog. It rises at every timestep 1 through 1000 (t=0.01 to 10). Initial factor at t=0 is 1. Final factor %.12f.\n',c(end));
fprintf(fid,'The energy plot shows the original saved samples, with red points every 100 steps for readability; these are not isolated contractility jumps. Full increment data are in the CSV.\n');
fprintf(fid,'Purse-string tension is reconstructed from the original recruitment recurrence, validated against saved final lambda and against the saved contractility trajectory. Final tension %.12f.\n',lambda(end));
fprintf(fid,'Perimeter energy coefficient = contractility factor / rstiff, so the bulk coefficient is 2 and final margin coefficient is %.12f. Physical time units have not been calibrated.\n',c(end)/s.param.rstiff);
fclose(fid);
copyfile('/tmp/plot_wound_contractility.m',fullfile(folder,'plot_contractility_annotations.m'));
fprintf('Increases: %d/%d; first step=%d; last step=%d; final contractility=%.12f; final lambda=%.12f\n',nnz(increased),n,find(increased,1),find(increased,1,'last'),c(end),lambda(end));
fprintf('Saved annotated energy, contractility plots, editable FIGs, per-step CSV and notes in %s\n',folder);
