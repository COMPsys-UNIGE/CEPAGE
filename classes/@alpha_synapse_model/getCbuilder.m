function str = getCbuilder(object)
    str = sprintf('alpha_synapse_model(%e,%e,%e,%e)',object.al,object.beta,object.nu,object.theta);
end