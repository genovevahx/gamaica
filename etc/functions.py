import datetime
from astropy.table import QTable
import matplotlib.pyplot as plt

def save_and_see(result):
    '''
    WARNING:
    astropy v3.2.1: The unit 'electron' can not be saved to FITS format
    therefore saving to ascii.ecsv for now

    ToDo:
    
    add ETC versioning information to meta data of table
    '''

    if not isinstance(result,QTable):
        raise TypeError('Input parameter must be a QTable instance')

    
    c=datetime.datetime.now()
    
    name = f"gamaica_etcrun_{c.hour:02.0f}h{c.minute:02.0f}m{c.second:02.0f}s_{c.year:04.0f}-{c.month:02.0f}-{c.day:02.0f}.dat"
    print('saving table to ',name)
    result.write(name,format='ascii.ecsv')

    ## visualize
    fig,axs=plt.subplots(2,2)
    fig.suptitle(result['exptime'][0])
    
    axs[0][0].plot(result['wavelength'],result['object'])
    axs[0][0].set_xlabel('AA')
    axs[0][0].set_ylabel(result['object'].unit)
    axs[0][0].set_title('Object flux')
    axs[0][0].semilogy(True)
    axs[0][0].set_xlim(result['wavelength'].value.min(),result['wavelength'].value.max())

    axs[0][1].plot(result['wavelength'],result['sky'])
    axs[0][1].set_xlabel('AA')
    axs[0][1].set_ylabel(result['sky'].unit)
    axs[0][1].set_title('Sky flux')
    axs[0][1].semilogy(True)
    axs[0][1].set_xlim(result['wavelength'].value.min(),result['wavelength'].value.max())


    axs[1][0].plot(result['wavelength'],result['noise'])
    axs[1][0].set_xlabel('AA')
    axs[1][0].set_ylabel(result['noise'].unit)
    axs[1][0].set_title('Total noise')
    axs[1][0].semilogy(True)
    axs[1][0].set_xlim(result['wavelength'].value.min(),result['wavelength'].value.max())

    axs[1][1].plot(result['wavelength'],result['snr'])
    axs[1][1].set_xlabel('AA')
    axs[1][1].set_ylabel(result['snr'].unit)
    axs[1][1].set_title('S/N')
    axs[1][1].semilogy(True)
    axs[1][1].set_xlim(result['wavelength'].value.min(),result['wavelength'].value.max())

    plt.subplots_adjust(left=0.12,bottom=0.11,right=0.99,top=0.93,wspace=0.26,hspace=0.42)

    figname = name.replace('.dat','.pdf')
    fig.savefig(figname,fmt='pdf',bbox_inches='tight')
    print('saving plots to ',figname)
    
    plt.show()

    
def view(tab):
    ''' This is not working '''
    if not isinstance(tab,QTable):
        raise TypeError('Input parameter must be a QTable instance')
    
    for c in tab.colnames:
        tab[c].format = '%.2f'
    tab.pprint(max_width=-1)
    

def list_templates():
    from etc.stddata import seddir
    '''List all available templates

    This function lists all templates that are available in the SED directory. 
    These are the only templates that can be used (for now)
    

    Returns
    -------
    list
        List of template names

    Examples
    --------

        >>> list_templates()[:3]
        ['Kinney_bulge', 'Kinney_ell', 'Kinney_s0']

    '''
    return sorted(map(lambda x: x.stem,seddir.glob('*/*.fits')))

