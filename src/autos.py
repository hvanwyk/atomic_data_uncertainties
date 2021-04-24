import re, numpy, pyatomdb, os
"""
Package for reading in AUTOSTRUCTURE files

Developed as required

usage:

a = autos.read_oic(dirname)

dirname is the directory where the oic file is. Should just be called oic.


Returns a dict of information from the file.

first index is NV
second index in LV
then there are entries for levels (lev)
then there are entries for configuration (cfg)
then there are entries for autoionization transition (aut)
then there are entries for radiative transition (trn)

which are all arrays which hopefully make sense.

This is easy to develop, I have only tested on a few files so there will be bugs.
Some of the stuff populating the "setup" doesn't make sense to me any more but
I'm sure I had my reasons.


"""


def increment_next_letters(let):
  if let[1] == 'z':
    let = chr(ord(let[0])+1)+'a'
  else:
    let = let[0]+chr(ord(let[1])+1)
  return let

def read_oic(directory):
  """
  Read the oic file. Directory is *either* the filename or the directory with the oic
  file in it.

  Parameters
  ----------
  directory:str
    Where the oic file lies

  Returns
  -------
  Structure with the data in it!

  """

  if os.path.isfile(directory)
    oicfile = directory
  else:
    oicfile = "%s/oic"%(directory)
  mode = 'CFG1'
  nextletters = 'aa'

  aut = {}
  ret = {}


  for line in open(oicfile, 'r'):
    if mode=='CFG1':
      if 'NV' in line:
        l = line.split()
        aut['setup']={}
        aut['setup']['NV']=int(l[1])
        aut['setup']['LV']=int(l[3])
        nshl = len(l)-5
        aut['K']= numpy.zeros(nshl, dtype=numpy.dtype({'names':['K','N','L'],\
                                                       'formats':[int, int, int]}))
      elif 'NL' in line:
        l = line.split()
        ncfg = int(line[0:3])
        Z= int(line[15:17])
        nel = int(line[23:26])
        l = line[30:].split()
        for il in range(len(l)):
          #print il, l[il]
          if il/2 >=nshl:
            aut['K']= numpy.append(aut['K'], \
                                 numpy.zeros(nshl, \
                                 dtype=numpy.dtype({'names':['K','N','L'],\
                                                    'formats':[int, int, int]})))
          if il %2==0:
            aut['K']['K'][il//2]=1+il/2
            aut['K']['N'][il//2]=int(l[il])
          else:
            aut['K']['L'][il//2]=int(l[il])
        aut['K']=aut['K'][:numpy.max([nshl, il//2])]
        mode = 'CFG2'
        aut['cfg']={}
        aut['cfg']=numpy.zeros(ncfg, dtype={'names':['CF','G','A','B','EIS',\
                                                     'cfgstr', 'parity'],\
                                            'formats':[int, int, int, \
                                                       int, '|30S', '|S40',\
                                                       int]})
        icfg = 0
    elif mode=='CFG2':
      if icfg < ncfg:
        l = line.split()
        aut['cfg']['CF'][icfg]=int(line[:5])
        aut['cfg']['G'][icfg]=int(line[5:10])
        aut['cfg']['A'][icfg]=int(line[13:15])
        aut['cfg']['B'][icfg]=int(line[15:17])

        aut['cfg']['EIS'][icfg]=line[19:].strip()
        #print(aut['cfg']['EIS'][icfg])
        if "*" in aut['cfg']['EIS'][icfg].decode():
          print("FIXING %s with %s"%(aut['cfg']['EIS'][icfg].decode(), nextletters))
          aut['cfg']['EIS'][icfg]=re.sub('\*', nextletters, aut['cfg']['EIS'][icfg].decode())

          nextletters = increment_next_letters(nextletters)
        #print(aut['cfg']['EIS'][icfg])
        tmp = pyatomdb.atomic.parse_eissner(aut['cfg']['EIS'][icfg], nel=nel)
        print("Parsed to %s"%(tmp))
        tmp2, x = pyatomdb.atomic.config_to_occup(tmp,\
                                            nel=nel)
        aut['cfg']['cfgstr'][icfg]=pyatomdb.atomic.occup_to_cfg(tmp2)
        aut['cfg']['parity'][icfg]=\
              pyatomdb.atomic.get_parity(aut['cfg']['cfgstr'][icfg])
        icfg+=1
      else:
        if 'NLEVEL' in line:
          mode='LEV1'
          nlev =  int(line[10:18])
          aut['lev'] = numpy.zeros(nlev, dtype=\
                         numpy.dtype({'names':['k','lv','t','s2p1','l',\
                                              'j2','cf','cfgstr','e_ryd',\
                                              'energy', 'parity'],\
                                      'formats':[int, int, int, int, int, \
                                                 int, int, '|S40', float, \
                                                 float, int]}))
          ilev = 0
        elif 'AUTO-IONIZATION' in line:
          mode='AUT1'
          
          aut['aut'] = numpy.zeros(10000, dtype=\
                         numpy.dtype({'names':['cfi','lvi','wi','cfc','lvc',\
                                              'wc','aa','ec_ryd','ei_ryd'],\
                                      'formats':[int, int, int, int, int, \
                                                 int, float, float, float]}))
          iauto = 0
          
    elif mode=='AUT1':
      mode = 'AUT2'
    elif mode=='AUT2':
      if 'NLEVEL' in line:
        mode='LEV1'
        aut['aut']=aut['aut'][:iauto]
        nlev =  int(line[10:18])
        
        aut['lev'] = numpy.zeros(nlev, dtype=\
                       numpy.dtype({'names':['k','lv','t','s2p1','l',\
                                            'j2','cf','cfgstr','e_ryd',\
                                            'energy', 'parity'],\
                                    'formats':[int, int, int, int, int, \
                                               int, int, '|S40', float, \
                                               float, int]}))
        ilev = 0
        
        continue



      l = line.split()
      if iauto >= len(aut['aut']):
          aut['aut']=numpy.append(aut['aut'], numpy.zeros(10000, dtype=aut['aut'].dtype))

      if len(l)==1:
        aut['aut'][iauto]['cfi']= -9999
        aut['aut'][iauto]['lvi']= -9999
        aut['aut'][iauto]['wi'] = -9999
        aut['aut'][iauto]['cfc']= -9999
        aut['aut'][iauto]['lvc']= -9999
        aut['aut'][iauto]['wc']=-9999
        aut['aut'][iauto]['aa']=-9999
        aut['aut'][iauto]['ec_ryd']=-9999
        aut['aut'][iauto]['ei_ryd']=float(l[0])
        
        

      else:
        aut['aut'][iauto]['cfi'] = int(l[0])
        aut['aut'][iauto]['lvi'] = int(l[1])
        aut['aut'][iauto]['wi']  = int(l[2])
        aut['aut'][iauto]['cfc'] = int(l[3])
        aut['aut'][iauto]['lvc'] = int(l[4])
        try:
          aut['aut'][iauto]['wc'] = int(l[5])
        except:
          if l[5]=='X':
            aut['aut'][iauto]['wc']=-9999
          else:
            raise()
        aut['aut'][iauto]['aa']    = float(l[6])
        aut['aut'][iauto]['ec_ryd']= float(l[7])
        
        try:
          aut['aut'][iauto]['ei_ryd']= float(l[8])
        except:
          print(line, l)
          raise()
          
      iauto+=1
      
      
            
    elif mode=='LEV1':

      if ilev <nlev:
        if 'RY' in line: continue
        l = line.split()
        aut['lev']['k'][ilev] = int(line[:5])
        aut['lev']['lv'][ilev] = int(line[5:10])
        aut['lev']['t'][ilev] = int(line[10:15])
        aut['lev']['s2p1'][ilev] = abs(int(line[15:20]))
        aut['lev']['l'][ilev] = int(line[20:25])
        aut['lev']['j2'][ilev] = int(line[25:30])
        aut['lev']['cf'][ilev] = int(line[30:35])
        aut['lev']['e_ryd'][ilev] = float(line[35:])
        icfg = numpy.where(aut['cfg']['CF']==aut['lev']['cf'][ilev])[0]
        aut['lev']['cfgstr'][ilev] = aut['cfg']['cfgstr'][icfg[0]]
        aut['lev']['energy'][ilev] = 1000*pyatomdb.const.RYDBERG*aut['lev']['e_ryd'][ilev]
        aut['lev']['parity'][ilev] = aut['cfg']['parity'][icfg[0]]
        ilev += 1
      else:
        if 'RADIATIVE' in line:
          mode='RAD1'
          aut['trn']=numpy.zeros(10000, dtype=numpy.dtype(
                                              {'names':['upper_cf','upper_lev',\
                                                        'upper_w','lower_cf',\
                                                        'lower_lev','lower_w',\
                                                        'ar','de_ryd'],\
                                               'formats':[int, int,\
                                                          int, int,\
                                                          int, int,\
                                                          float, float]}))
          # get the Z and N
          Zind = line.index('Z=')+2
          Z=int(line[Zind:Zind+2])
          Nind = line.index('N=')+2
          N=int(line[Nind:Nind+2])


          linespl = line.split()
          #try:
          #  for llll in linespl:
          #    if llll[0:2] == 'Z=':
          #      Z = int(llll[2:])
          #    if llll[0:2] == 'N=':
          #      N = int(llll[2:])
          #except:
          #  print(line)
          #  print(llll)
          #  raise()


          aut['setup']['Nelectrons'] = N
          aut['setup']['Z'] = Z


          itrn = 0

    elif mode=='RAD1':
      if 'DEL' in line: continue
      if len(line) < 10:
        aut['trn']=aut['trn'][:itrn]
        
        if not aut['setup']['NV'] in ret.keys():
          ret[aut['setup']['NV']]={}

        ret[aut['setup']['NV']][aut['setup']['LV']]=aut
        mode='CFG1'  
      else:
        l = line.split()
        if itrn >=len(aut['trn']):
          aut['trn'] = numpy.append(aut['trn'], numpy.zeros(10000, dtype=numpy.dtype(
                                              {'names':['upper_cf','upper_lev',\
                                                        'upper_w','lower_cf',\
                                                        'lower_lev','lower_w',\
                                                        'ar','de_ryd'],\
                                               'formats':[int, int,\
                                                          int, int,\
                                                          int, int,\
                                                          float, float]})))
        aut['trn']['upper_cf'][itrn]=int(line[:5])
        aut['trn']['upper_lev'][itrn]=int(line[5:10])
        aut['trn']['upper_w'][itrn]=int(line[10:15])
        aut['trn']['lower_cf'][itrn]=int(line[15:20])
        aut['trn']['lower_lev'][itrn]=int(line[20:25])
        aut['trn']['lower_w'][itrn]=int(line[25:30])
        aut['trn']['ar'][itrn]=abs(float(line[30:45]))
        aut['trn']['de_ryd'][itrn]=float(line[45:60])
        itrn+=1


  return ret





