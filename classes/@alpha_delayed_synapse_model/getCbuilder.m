function str = getCbuilder(object)
    str = sprintf('alpha_delayed_synapse_model(%e,%e,%e,%e,%e,%e)',object.al,object.beta,object.nu,object.theta,object.delays(1),object.g2);
end