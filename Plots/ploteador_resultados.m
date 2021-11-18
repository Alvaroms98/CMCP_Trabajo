%% Ploteo de los resultados del trabajo de CMCP 2.2: Poisson MPI-OpenMP

clear 
close all
clc

%% Carga de los datos
load('datos_trabajo_CMCP.mat')

%%%%%%%%%%%%%%%%%% OpenMP vs MPI %%%%%%%%%%%%%%%%%%

%% Speedup Fuerte 400x400
fig1 = figure(1);
fig1.Color = [1,1,1];

plot(OpenMP.Hilos, OpenMP.SpeedupFuerte400x400, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.SpeedupFuerte400x400, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.Procesos, 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('Strong scaling - 400x400')
ylabel('Speedup [-]')
xlabel('Processes/Threads')
legend('OpenMP','MPI','Ideal','Location','northwest')
xlim([1,64])
ylim([1,64])

%% Speedup Fuerte 1000x1000
fig2 = figure(2);
fig2.Color = [1,1,1];

plot(OpenMP.Hilos, OpenMP.SpeedupFuerte1000x1000, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.SpeedupFuerte1000x1000, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.Procesos, 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('Strong scaling - 1000x1000')
ylabel('Speedup [-]')
xlabel('Processes/Threads')
legend('OpenMP','MPI','Ideal','Location','northwest')
xlim([1,64])
ylim([1,64])


%% Speedup Débil

fig3 = figure(3);
fig3.Color = [1,1,1];

plot(OpenMP.Hilos, OpenMP.SpeedupDebil, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.SpeedupDebil, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.Procesos, 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('Weak scaling')
ylabel('Speedup [-]')
xlabel('Processes/Threads')
legend('OpenMP','MPI','Ideal','Location','northwest')
xlim([1,64])
ylim([1,64])

%% Eficiencia Fuerte 400x400

fig4 = figure(4);
fig4.Color = [1,1,1];

plot(OpenMP.Hilos, OpenMP.EficienciaFuerte400x400, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.EficienciaFuerte400x400, 'LineWidth', 3);
hold on
[tam, ~] = size(MPI.Procesos);
plot(MPI.Procesos, 100*ones(tam,1), 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('Strong scaling - 400x400')
ylabel('Efficiency [%]')
xlabel('Processes/Threads')
legend('OpenMP','MPI','Ideal','Location','southwest')
ylim([0,120])
xlim([1,64])

%% Eficiencia Fuerte 1000x1000

fig5 = figure(5);
fig5.Color = [1,1,1];

plot(OpenMP.Hilos, OpenMP.EficienciaFuerte1000x1000, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.EficienciaFuerte1000x1000, 'LineWidth', 3);
hold on
[tam, ~] = size(MPI.Procesos);
plot(MPI.Procesos, 100*ones(tam,1), 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('Strong scaling - 1000x1000')
ylabel('Efficiency [%]')
xlabel('Processes/Threads')
legend('OpenMP','MPI','Ideal','Location','southwest')
ylim([0,120])
xlim([1,64])

%% Eficiencia Débil

fig6 = figure(6);
fig6.Color = [1,1,1];

plot(OpenMP.Hilos, OpenMP.EficienciaDebil, 'LineWidth', 3);
hold on
plot(MPI.Procesos, MPI.EficienciaDebil, 'LineWidth', 3);
hold on
[tam, null] = size(MPI.Procesos);
plot(MPI.Procesos, 100*ones(tam,1), 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('Weak scaling')
ylabel('Efficiency [%]')
xlabel('Processes/Threads')
legend('OpenMP','MPI','Ideal','Location','southwest')
ylim([0,120])
xlim([1,64])


%%%%%%%%%%%%%%%%%% MPI-OpenMP %%%%%%%%%%%%%%%%%%

%% Speedup

fig7 = figure(7);
fig7.Color = [1,1,1];

plot(OpenMP_MPI.Procesadores, OpenMP_MPI.SpeedupFuerte, 'LineWidth', 3);
hold on
plot(OpenMP_MPI.Procesadores, OpenMP_MPI.SpeedupDebil, 'LineWidth', 3);
hold on
plot(OpenMP_MPI.Procesadores, OpenMP_MPI.Procesadores, 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('MPI-OpenMP')
ylabel('Speedup [-]')
xlabel('Processes/Threads')
legend('Fuerte 800x800','Debil','Ideal','Location','northwest')
xlim([1,256])
ylim([1,20])


%% Eficiencia 

fig8 = figure(8);
fig8.Color = [1,1,1];

plot(OpenMP_MPI.Procesadores, OpenMP_MPI.EficienciaFuerte, 'LineWidth', 3);
hold on
plot(OpenMP_MPI.Procesadores, OpenMP_MPI.EficienciaDebil, 'LineWidth', 3);
hold on
[tam,~] = size(OpenMP_MPI.Procesadores);
plot(OpenMP_MPI.Procesadores, 100*ones(tam,1), 'Color', [0,0,0], 'LineWidth', 1);

grid on
title('MPI-OpenMP')
ylabel('Efficiency [%]')
xlabel('Processes/Threads')
legend('Fuerte 800x800','Debil','Ideal','Location','northwest')
xlim([1,256])
ylim([0,120])


%%%%%%%%%%%%%%%%% Guardar figuras %%%%%%%%%%%%%%%%%%%%%%%%%%

%% Guardar imágenes

saveas(fig1, "speedup_fuerte_400x400.png");
saveas(fig2, "speedup_fuerte_1000x1000.png");
saveas(fig3, "speedup_debil_open_vs_mpi.png");
saveas(fig4, "eficiencia_fuerte_400x400.png");
saveas(fig5, "eficiencia_fuerte_1000x1000.png");
saveas(fig6, "eficiencia_debil_open_vs_mpi.png");
saveas(fig7, "speedup_MPI_y_OpenMP.png");
saveas(fig8, "eficiencia_MPI_y_OpenMP.png");
    