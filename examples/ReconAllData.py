import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from runReconstruction import runReconstruction

configfiles = ['config/config_kidneytissue.yml', 'config/config_u2os_cells.yml']
#TODO: test analyzing sub-z stack on Kidney Tissue data.

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

