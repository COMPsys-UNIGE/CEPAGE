function str = getCbuilder(object)
    str = sprintf('A_synapse_model(%e,%e,%e,%e)',object.al,object.beta,object.nu,object.theta);
end