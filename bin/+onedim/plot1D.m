function [h] = plot1D(mesh,f)

z = mesh.z;
h = plot(z,f,'.-');
xlim([min(mesh.z),max(mesh.z)])
ylim([-2,2])
hold on

end

