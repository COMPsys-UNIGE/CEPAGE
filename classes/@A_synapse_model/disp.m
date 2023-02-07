function disp(object)

disp('Andrey synapse model.');
disp('The synapse activation function is');
disp('   ');
disp('x');
disp('   ');
disp('The state equations is:');
disp(' alpha*s*(1-s)*1/(1+exp(-nu*(Vpre-theta))-beta*s');

disp(' ');
disp('Parameters value:');
disp(['     alpha = ',num2str(object.al)]);
disp(['      beta = ',num2str(object.beta)]);
disp(['     theta = ',num2str(object.theta)]);
disp(['     nu    = ',num2str(object.nu)]);

end