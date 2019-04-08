from workflow.runReconstruction import runReconstruction

configfiles = ['config/config_u2os_cells.yml','config/config_kidneytissue.yml']
#TODO: analysis runs fine on glass beads and u2os cells, but not on Kidney Tissue data and Mouse Brain slice. Needs debugging.

if __name__ == '__main__':
    """
    Reconstruct all data supplied with ReconstructOrder repository.
    The total size of dataset is approximately 4GB and will be available for downloaded from the release on github.
    The reconstruction parameters are specified in the configuration files stored in config folder.

    Parameters
    ----------
    Path of config files relative to the root of the repository.
    
    Returns
    --------
    Outputs data to disk.
    """

    for configfile in configfiles:
        print(configfile + '\n--------------')
        runReconstruction(configfile)
        print('\n--------------')


