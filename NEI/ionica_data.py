'''
This is a data storage location for Ionica by Marcus Dupont

Atomic Data used: AtomDB by: Adam Foster

'''



def read_atomdb_data(elements=['H', 'He', 'C',     # twelve most abundant elements
                               'N', 'O', 'Ne',
                               'Mg', 'Si', 'S', 
                               'Ar', 'Ca', 'Fe', ] , 
                     data_directory= '/home/mdupont/atomdb/APED', 
                     screen_output=False):
    return elements