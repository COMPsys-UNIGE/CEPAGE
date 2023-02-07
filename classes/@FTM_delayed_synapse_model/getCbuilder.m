function str = getCbuilder(object)
    str = sprintf('FTM_delayed_synapse_model(%e,%e,%e,%e,%e)',object.nu,object.theta,object.g1,object.g2,object.delays(1));
end