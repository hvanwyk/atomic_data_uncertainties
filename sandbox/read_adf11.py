# python3 function to read adf11 file
# returns numpy array with shape (z_max, ntemps, ndens)
# note that array is zero indexed, so z_max=1 is index 0
# Author: Ivan Arnold

def read_adf11(filename):
    import numpy as np  
    with open(filename,'r') as f:  # open the file
        
        line_dat = f.readline()  # read first line
        param = line_dat.split()
        z_max = int(param[0])    # max ion charge
        ndens = int(param[1])    # number of densities
        ntemps = int(param[2])   # number of temperatures
        
        # determine the number lines in each density block
        dens_quot, dens_rem =  divmod(ndens,8)
        dens_lines = dens_quot
        if dens_rem > 0:
          dens_lines = dens_lines +1
        
        #determine the number of lines in each temperature block
        temp_quot, temp_rem =  divmod(ntemps,8)
        temp_lines = temp_quot
        if temp_rem > 0:
          temp_lines = temp_lines + 1
          
        
        
        while line_dat[1] != '-':  # skip past first line that begins with '-'
          line_dat = f.readline()
        
        # create an array with the density values
        dens_array = np.array(0)
        for i in range(dens_lines):
            dat = str.strip(f.readline()).split()
            dat = np.array(dat)
            dens_array = np.hstack((dens_array,dat))
            dens_array =  dens_array.astype(float)
            dens = 10 ** (dens_array[1:]) # standard ADAS format uses log10
        
        # create an array with the temperature values
        temp_array = np.array(0)
        for i in range(temp_lines):
            dat = str.strip(f.readline()).split()
            dat = np.array(dat)
            temp_array = np.hstack((temp_array,dat))
            temp_array =  temp_array.astype(float)
            temps = 10 ** (temp_array[1:]) # convert from log10
        
        # create a data array with shape ((z_max, ntemps, ndens))
        # and populate with the GCR coefficients from adf11 file
        adf11_dat = np.zeros((z_max, ntemps, ndens))
        for k in range(z_max):
            f.readline()
            qcd_array = dens_array
            for i in range(ntemps):
              ldat = qcd_array
              cdat = np.array(0)
              for j in range(dens_lines):
                dat = str.strip(f.readline()).replace('D','E').split()
                dat = np.array(dat)
                cdat = np.hstack((cdat,dat))
              qcd_array = np.vstack((ldat,cdat))
              qcd = (qcd_array.astype(float)[1:,1:])
            adf11_dat[k, :, :] = 10 ** qcd
        return(dens, temps, adf11_dat)


if __name__ == "__main__":
    dens, temps, acd_dat = read_adf11('acd96_he.dat')
    dens, temps, scd_dat = read_adf11('scd96_he.dat')
    print(dens)
                  
                  

