import subprocess
import os
import threading
#import pandas as pd
#import matplotlib.pyplot as plt

class epipy:
    def __init__(self,solute_structure,solute_topology,to_gro=False,gen_idc=False,
                 convert=False,recenter=False,box_size=[10,10,10]):
      """
      EPIPY: Mk.II.beta
      3DRISM interface ,water placement and data manipulation tool
      =============================================
      class is initialized with a solute structure file and its corresponding topology
      +++++++++++++++++++++++++++++++++++++++++++++
      solute_structure :: .gro, .pdb or .xtc file
      solute_topology :: .top or .solte file
      gen_idc :: bool, generate an ion-dipole-correction (IDC) enabled .solute file
      this file will be used for the remainder of the calculation
      convert :: convert solute_topology to .solute file

      ========== IF GROMACS ====================
      - requirements: GROMACS installed and properlly sourced
      to_gro :: convert solute_structure file to .gro using pdb2gmx
      recenter :: recenter solute incenter of box using gmx -edifconf

      """
      self.cmd_path = "" #path to eprism command
      self.grid = [] # number of grids per unit cell LxWxH
      self.closure = 'PSE3'
      self.path = ""
      self.ndiis = 5
      self.solute_path = ""
      self.solvent_path = ""#~/EPISOL-1.1.326/solvent/"
      self.r_c = 1
      self.coulomb = 'coulomb'
      self.delvv = 0.5
      self.err_tol = 1e-08
      self.dynamic_delvv = 1
      self.T = 298
      self.log = ''#'episol_out.log'
      self.threads = threading.active_count()
      self.rism_args = str()
      self.test_args = str()
      self.rism_cmd = str()
      #self.out_file = self.log#str()
      self.get_eprism_path = ''.join([chr(i) for i in subprocess.check_output(['whereis eprism3d'],shell=True)[:-1]]).split()[1][:-8]
      ## initialize solute and solvent right at the begining
      # you can do this later onif you want to change files
      self.solute(solute_structure,solute_topology,to_gro,gen_idc,convert,recenter,box_size)
      self.solvent()
      # here we set the default report name to the solute structure file name
      # and the default save to output ALL values in the log, there is no reason
      # not to do this as it is like 2 lines of text
      self.rism(step=500,resolution=1,args=('all'))
      self.report(out_file_name=f'{self.structure_file[:-4]}_out',args=('all'))
      ####################################################################################
      # HERE I AM MOVING THIS TO __INIT__
      ####################################################################################
    def solute(self,solute_structure,solute_topology,to_gro=False,gen_idc=False,convert=False,recenter=False,box_size=[10,10,10]):
        """
        solute_structure :: .gro, .pdb or .xtc file
        solute_topology :: .top or .solte file
        gen_idc :: bool, generate an ion-dipole-correction (IDC) enabled .solute file
        this file will be used for the remainder of the calculation
        convert :: convert solute_topology to .solute file

        ========== IF GROMACS ====================
        - requirements: GROMACS installed and properlly sourced
        to_gro :: convert solute_structure file to .gro using pdb2gmx
        recenter :: recenter solute incenter of box using gmx -edifconf
        """

        # I wanted to get the file path without using the os package as it does NOT
        # I REPEAT NOT work well for distributed computing so please dont hate me for
        # what youre about to see
        self.structure_file = solute_structure.split('/')[-1] # get rid of the path
        if solute_structure.find('/') != -1:
            self.solute_path = '/'.join(solute_structure.split('/')[:-1])+"/" # this is goofy
        else:
            self.solute_path = ''

        self.solute_top = solute_topology.split('/')[-1]
        if solute_topology.find('/') != -1:
            self.solute_top_path = '/'.join(solute_topology.split('/')[:-1])+"/"
        else:
            self.solute_top_path = ''
        # resume
        if gen_idc:
            convert = True

        if convert:
            subprocess.run([f"gmxtop2solute -p {self.solute_path+self.solute_top} -o {self.solute_top[:-4]}.solute"],shell=True)
            print(f"converted {self.solute_top} to {self.solute_top[:-4]}.solute")
            self.solute_top = f"{self.solute_top[:-4]}.solute"
            """ We will reset our path since we just wrote a new fle to the current dir"""
            self.solute_top_path = ''
        if gen_idc:
            #subprocess.run([f"{self.get_eprism_path}generate-idc.sh {self.solute_top[:-4]}.solute > idc_{self.solute_top[:-4]}.solute"],shell=True)
            subprocess.run([f"generate-idc.sh {self.solute_top} > idc_{self.solute_top}"],shell=True)
            self.solute_top = f'idc_{self.solute_top}'
            print(f"generated idc-enabled solute file to: {self.solute_top}")
        #else:
        #    self.solute_top = solute_topology.split('/')[-1]
        ####################################################################################
        # GROMACS REQUIRED FOR THESE COMMANDS, may remove later on but can be handy
        ####################################################################################
        if recenter:
            subprocess.run([f"gmx  editconf -f {self.structure_file} -c yes -box {box_size[0]} {box_size[1]} {box_size[2]} -o {self.structure_file[:-4]}.gro"],shell=True)
            self.structure_file = f"{self.structure_file[:-4]}.gro"
        if to_gro:
            subprocess.run([f"gmx  editconf -f {self.structure_file} -o {self.structure_file[:-4]}.gro"],shell=True)
            self.structure_file = f"{self.structure_file[:-4]}.gro"
        #box = [box_size,box_size,box_size]
        #self.rism_args += f" -f {self.path+self.solute_path+self.structure_file} -s {self.path+self.solute_path+self.solute_top}"
      ####################################################################################
      # this block reads the structure file and returns the box dimensions
        if f"{self.path+self.solute_path+self.structure_file}"[-3:] == 'pdb':
          self.file_type = 'pdb'
          with open(f"{self.path+self.solute_path+self.structure_file}",'r') as sol:
              for line in sol:
                tmp = line.split()
                if tmp[0] == "CRYST1":
                  self.solute_box = [tmp[1],tmp[2],tmp[3]]
                  break
          sol.close()

        elif f"{self.path+self.solute_path+self.structure_file}"[-3:] == 'gro':
          self.file_type = 'gro'
          """HELP FROM:
          https://stackoverflow.com/questions/3346430/
          what-is-the-most-efficient-way-to-get-first-and-last-line-of-a-text-file/18603065#18603065"""
          import os
          with open(f"{self.path+self.solute_path+self.structure_file}",'rb') as f:
              try:
                  f.seek(-2, os.SEEK_END)
                  while f.read(1) != b'\n':
                      f.seek(-2, os.SEEK_CUR)
              except OSError:
                  f.seek(0)
              last_line = f.readline().decode()
          f.close()
        self.solute_box = [float(val) for val in last_line.split()]

    def solvent(self,solvent_topology=None):
        """
        Reads the solvent topology file.
        If no file specified, will default search
        the site-packages dir where episol.epipy is
        stored. this may fail as the os package sometimes
        encounters errors
        =================================
        solvent_topology :: solvent correlation file
        -------------------------------
        if solvent topology is set to None i.e. has no input
        then we will select the tip3p 0.01A file from the source directory
        of EPIPY, this way to avoid annoying paths
        """
        from os.path import dirname
        from inspect import getfile
        # here we will set the solvent topology default to search the
        # site-packages directory. this way not as troublesome
        if not solvent_topology:
         self.solvent_top ='tip3p-amber14.01A.gaff'
         self.solvent_path = f'{dirname(getfile(epipy))}/'
        else:
            self.solvent_top = solvent_topology
        #self.rism_args += f" -p {self.path+self.solvent_path+self.solvent_top}"
        ####################################################################################
        ####################################################################################
    def rism(self,step=500,resolution=1,args=('all')):
        """
        Sets the steps and resolution of the calculation.
        Will automatically set number of grids based on the box dimensions
        acquired from self.solute()
        =======================================================
        step :: int(), Number of SCF steps to perform
        resolution :: int(), grid resolution of box (will override if previously set)
        args :: string(), values to save to the dump file; options below:
        ++++++++++++++++++++++++++ SAVE ARGS ++++++++++++++++++++++++++
        all : reports all of the below
        command : save the exact command that was run when the file was produced
        guv : g(r) foreach grid point
        ld : smoothed g(r)  i.e. g(r) convolved with kernel
        coul : coulombic potential at each grid
        excess : mu^ex at each grid
        """
        self.rism_step = step
        self.resolution = resolution
        #self.rism_cmd +=  f" -cmd closure={self.closure} rism,step={self.rism_step}"#f" rism,step={self.rism_step}"
        ###
        tmp = self.solute_box
        if self.file_type == 'gro':
          # box is in nm
          convert = (10)/resolution
          self.grid = [i*convert for i in tmp]
        if self.file_type == 'pdb':
          # box is in Angstrom
          convert = 1/resolution
          self.grid = [i*convert for i in tmp]

        self.save_command = f' save:{args}'

    def report(self,out_file_name:str,args=('all')):
        """
        function specifies the name of the .log file and
        which params to write out into it
        file name can be specified with self.log as well
        ================================================
        out_file :: name for .log file
        **args :: command save strings
        """
        self.log = out_file_name
        cmd_string = " report:"
        for arg in args:
            cmd_string += f"{arg}"
        self.to_report = f" {cmd_string}" #{self.rism_cmd} {cmd_string}"

    def get_version(self):
        """
        returns the version by simply calling episol kernel
        from the command line
        """
        temp = subprocess.check_output([f'{self.cmd_path}eprism3d --version'],shell=True)
        return ''.join([chr(i) for i in temp[:-1]])

    def get_help(self,search_str:str):
      """
      searches episol by calling --help from the CLI
      ==============================================
      search_str :: string to search for
      """
      search_str = search_str.strip()
      xx = subprocess.run(["eprism3d","--h",search_str],capture_output=True)
      xv = ''.join([chr(i) for i in xx.stdout[:-1]])
      print(xv)

    def test(self,nt=1,v=1):
        """
        Calls episol and adds -test flag to the CLI
        ==========================================
        v :: int(), verbose rating
        nt :: int(), number of threads to use
        """
        self.test_args += f" -f {self.path+self.solute_path+self.structure_file} -s {self.path+self.solute_path+self.solute_top}"
        self.test_args += f" -p {self.path+self.solvent_path+self.solvent_top}"
        self.test_args += f" -coulomb {self.coulomb}"
        #self.rism_cmd += f' -cmd closure={self.closure} rism,step={self.rism_step}'
        self.test_args += f" -rc {self.r_c} "
        self.test_args += f" -T {self.T}"
        self.test_args += f" -log {self.log}.log"
        self.test_args += f" -nr {self.grid[0]}x{self.grid[1]}x{self.grid[2]}"
        self.test_args += f" -ndiis {self.ndiis}"
        self.test_args += f" -errtolrism {self.err_tol}"
        self.test_args += f" -delvv {self.delvv}"
        self.test_args += f" -dynamic-delvv {self.dynamic_delvv}"
        self.test_args += " -pwd ./" #"~/mnt/f/water_proj"
        self.test_args += f" -o {self.log}"
        self.test_args += f" -cmd closure={self.closure} rism,step={self.rism_step}"#f" rism,step={self.rism_step}"#self.rism_cmd
        self.test_args += self.to_report
        self.test_args += self.save_command
        ###
        self.test_args += f' -nt {nt} -test'
        subprocess.run([f"{self.cmd_path}eprism3d {self.test_args}"],shell=True)
        self.test_args = ''

        with open(f"{self.log}.log",'r') as f:
          #mem = str()
          for line in f:
            tmp = line.split()
            if tmp[0] == "Memory":
              mem = f'{tmp[5]}{tmp[6]}'
              break
        f.close()
        return mem

    def kernel(self,nt=1,v=1):
        """
        v :: int(), verbose rating
        nt :: int(), number of threads to use
        """
        self.rism_args += f" -f {self.path+self.solute_path+self.structure_file} -s {self.path+self.solute_top_path+self.solute_top}"
        self.rism_args += f" -p {self.path+self.solvent_path+self.solvent_top}"
        self.rism_args += f" -coulomb {self.coulomb}"
        #self.rism_cmd += f' -cmd closure={self.closure} rism,step={self.rism_step}'
        self.rism_args += f" -rc {self.r_c} "
        self.rism_args += f" -T {self.T}"
        self.rism_args += f" -log {self.log}.log"
        self.rism_args += f" -nr {self.grid[0]}x{self.grid[1]}x{self.grid[2]}"
        self.rism_args += f" -ndiis {self.ndiis}"
        self.rism_args += f" -errtolrism {self.err_tol}"
        self.rism_args += f" -delvv {self.delvv}"
        self.rism_args += f" -dynamic-delvv {self.dynamic_delvv}"
        self.rism_args += " -pwd ./" #"~/mnt/f/water_proj"
        self.rism_args += f" -o {self.log}"
        self.rism_args += f" -cmd closure={self.closure} rism,step={self.rism_step}"#f" rism,step={self.rism_step}"#self.rism_cmd
        self.rism_args += self.to_report
        self.rism_args += self.save_command
        # using unix system is better to run subprocess as
        # argument strings, rather than **args
        # i.e. better "exe.exe -command value" than "exe.exe","command","value"
        #self.rism_args += f' save:all -nt {nt} -v {v}'
        self.rism_args += f' -nt {nt} -v {v}'
        subprocess.run([f"{self.cmd_path}eprism3d {self.rism_args}"],shell=True)#,self.rism_args],shell=True)
        self.rism_args = '' # clear vars so we do not overwrite if reinitializing
        with open(f'{self.log}.log','r') as f:
          for line in f:
            tmp = line.split()
            if tmp[0][:4] == "RISM":
              out = float(tmp[4])
              step = tmp[2]
        f.close()
        try:
            if out > self.err_tol:
                print(f"Failed to reach desired err_tol of {self.err_tol}")
                print(f"Actual error: {out}")
                print(f"Difference: {self.err_tol - out}")
                print(f"RISM finished at step {step}")
            else:
                print(f"Calculation finished in {step} steps ")
                print(f"err_tol: {self.err_tol} actual: {out} ")
        except UnboundLocalError:
            with open(f'{self.log}.log','r') as rr:
                tmp_err_hold = rr.read()
            rr.close
            print("It appears your calculation has encountered an error")
            print(f"Please see the output from the log file ({self.log}.log)")
            print(tmp_err_hold)
        return

    def dump(self,file_name='',out_name=False,list_values=False,value_to_extract=1):
        """
        Extracts the compressed .ts4s file and writes to a txt file
        ==================================
        file_name :: dump file to read
        out_name :: extract files to this txt file
        list_values :: return list of saved values
        value_to_extract :: given the list of values select index to extract
        """
        if not out_name:
            out_name = file_name

        if file_name == '':
          file_name = f'{self.log}.ts4s'

        if list_values:
            thold = ''.join([chr(i) for i in subprocess.check_output([f"ts4sdump -f {file_name} -l"],shell=True)]).split('\n')
            for item in thold:
              print(item)
        else:
          subprocess.run([f"ts4sdump -f {file_name} -e {value_to_extract} > {out_name}.txt"],shell=True)
        self.extracted_file = f'{out_name}.txt'

    def err(self,log_file_name=None):
      """
      This function reads a log file and returns
      the SCF stdev in an array
      ========================
      log_file_name :: string, name of .log file
      ======================== Returns
      out_arr :: np.array() with shape (1, # steps)
      array index represents the SCF step
      """
      from numpy import array

      if not log_file_name:
          log_file_name = f'{self.log}'
      with open(f'{log_file_name}.log','r') as f:
          out_arr = []
          for line in f:
              tmp = line.split()
              if tmp[0][:4] == "RISM":
                  out_arr.append(float(tmp[4]))
                  #step = tmp[2]
      f.close()
      #print(out_arr)
      return array(out_arr)

    def reader(self,file_in:str,laplacian:bool=False,convolve:bool=False,sigma:float=1.52,file_out:str='out',dx=False):
      """
      This function takes in an uncompressed dump file txt file
      it reads the grid size and shape
      If specified we can save/export to a dx file which can be loaded into pymol/vmd/etc.
      =================================
      file_in :: the decompress txt file from our ts4s dump command
      grid_spacing :: [x,y,z] values for \delta grid, e.g. grid spacing of 0.5A would
      have grid_spacing = [0.5,0.5,0.5]
      ---------- IF DX=TRUE ---------------
      file_out :: filename for saved dx file
      ################# WARNINGS ################
      since the dx file was made by IBM in the 90s
      many nont-so-modern softwares will struggle to read
      dx files with comments and it appears many follow their own
      format specifications
      """
      from numpy import loadtxt,zeros,copy
      from scipy.ndimage import laplace,gaussian_laplace

      grid_spacing = [self.resolution for _ in range(3)]
      #[i/j for (i,j) in zip(self.solute_box,self.grid)]

      xs,ys,zs = loadtxt(file_in)[-1][:3]
      x = loadtxt(file_in,usecols=(3))
      #xs,ys,zs
      shaped = zeros((int(xs),int(ys),int(zs)))
      cont = int()
      for zval in range(int(zs)):
        for yval in range(int(ys)):
            for xval in range(int(xs)):
                shaped[xval][yval][zval] = x[cont]#np.linalg.norm(x[cont])
                cont +=1
      if laplacian and not dx:
        return laplace(shaped)
      elif convolve and not dx:
          return gaussian_laplace(shaped,sigma=sigma)
      elif not dx:
        return shaped
      if dx:
        with open(f'{file_out}.dx','w+') as outfile:
          outfile.write(f"""object 1 class gridpositions counts     {xs}     {ys}      {zs}
origin 0.0000   0.0000   0.0000
delta  {grid_spacing[0]} 0 0
delta  0 {grid_spacing[1]} 0
delta  0 0  {grid_spacing[2]}
object 2 class gridconnections counts     {xs}     {ys}      {zs}
object 3 class array type double rank 0 items {int(xs*ys*zs)} follows\n""")
          for ind, val in enumerate(shaped.flatten()):
            outfile.write(f"{val:0.5e} ")
            if (ind != 0) and (ind % 3 == 0):
              outfile.write("\n")
        outfile.close()

    def placement(self,num_waters_to_place:int,radius=1.9269073728633292,filename=False,grid_spacing=None,write_pdb=False,outname='out.pdb',weight=None):
      import numpy as np
      import datetime
      """function selects top distribution
      then places water there, then removes that density
      and continues to look. kind of like placevent.
      =============================
      filename:: txt file containing density
      inputgrid:: md.Grid object
      outname:: string for pdb file
      num_waters_to_place:: number of waters to place
      radius:: distance around selected point to omit from remaining placement
      ideal_radius = ((10**3)/(4/3)/(np.pi)/33.3679)**(1/3) based on number density of water
      grid_spacing:: conversion from array indices to angstrom, units of 1/A
      """
      if grid_spacing is None:
        # you can overide the grid spacing if need be
        grid_spacing = self.resolution
      if not filename:
        filename = f'{self.log}.txt'

      xs,ys,zs = np.loadtxt(filename)[-1][:3]
      x = np.loadtxt(filename,usecols=(3))

      shaped = np.zeros((int(xs),int(ys),int(zs)))
      cont = int()
      for zval in range(int(zs)):
        for yval in range(int(ys)):
            for xval in range(int(xs)):
                shaped[xval][yval][zval] = x[cont]#np.linalg.norm(x[cont])
                cont +=1

      inputgrid = np.copy(shaped)
      x = np.arange(0,inputgrid.shape[0])
      y = np.arange(0,inputgrid.shape[1])
      z = np.arange(0,inputgrid.shape[2])
      out_array = []
      for wat in range(num_waters_to_place):
          r_x,r_y,r_z = np.where(inputgrid == np.max(inputgrid))
          r_x,r_y,r_z = r_x[0],r_y[0],r_z[0]
          out_array.append([float(r_x),float(r_y),float(r_z)])
          mask = (x[:,np.newaxis,np.newaxis]-r_x)**2 + (y[np.newaxis,:,np.newaxis]-r_y)**2 + (z[np.newaxis,np.newaxis,:]-r_z)**2  < (radius/grid_spacing)**2
          inputgrid[mask] = -1
      ################################
      if write_pdb:
        x,y,z = self.grid
        with open(f'{outname}.pdb','w+') as pdb:
            pdb.write(f"CRYST1    {str(x).ljust(4,'0')}    {str(y).ljust(4,'0')}    {str(z).ljust(4,'0')}  90.00  90.00  90.00 P 1           1\n")
            count = 1
            for val in np.array(out_array)*grid_spacing:#guess2:
                pdb.write(f"ATOM   {str(count).rjust(4,' ')}  O    HOH  A{str(count).rjust(3,' ')} {val[0]:8.3f} {val[1]:8.3f} {val[2]:8.3f}\n")
                count +=1

      return np.array(out_array)*grid_spacing

    def select_coords(self,in_file:str,sele:str,atom_sele:str=None,conv_fact:int=10):
      """
      in_file: .gro file
      sele: selection string -> resname only so far
      returns:

      ========
      conv_fact: convert from nm to grid-bits
      for .gro this is X nm*(10A/1nm)
      """
      import numpy as np
      # this is where we select the selection string
      tmp_string = sele.split()
       # Xnm*(10A/1nm)
      with open(f'{in_file}','r') as r:
          count = int()
          atom_count = int()
          out_dict = {}
          for line in r:
            count += 1
            if count == 2:
              atom_count = int(line.split()[0])
            if count == atom_count+2:
              break # this is the end of the file
            ###########
            if count > 2: # ignore header
              try:
                res_id = line[:5].split()[0]
                res_name = line[5:10].split()[0]
                # there is probably a better way of doing this
                if f"{res_name+res_id}" not in out_dict.keys():
                  # if the residue is not in the dictionary, add it
                  out_dict[f"{res_name+res_id}"] = {}
                else:
                  #print(line[20:28].split()[0])
                  out_dict[f"{res_name+res_id}"][f"{line[10:15].split()}"] = conv_fact*np.array([float(line[20:28].split()[0]),
                                                                                      float(line[28:36].split()[0]),
                                                                                      float(line[36:44].split()[0])])
              except ValueError:
                print(f"ERROR OCCURED AT LINE: {count}")
                continue
      # now we return the x,y,z positions ONLY
      out_array = np.array([])
      for val in out_dict.keys():
        if val[:3] == sele:
          # this uh is a pretty bad way of doing this
          out_array = np.append(out_array,[i for i in out_dict[val].values()])
      #out_array = np.array(out_array)
      #print('out array',out_array)
      """if len(out_array.shape) > 2:
        # this is a bad way of doing this
        y,z,_ = out_array.shape
        out_array = out_array.flatten()#.shape
        out_array = out_array.reshape((y*z,3))"""
      out_array = out_array.reshape((int(len(out_array)/3),3))

      return out_array

    def select_around(self,rism_grid,in_coords,around=5.0):
      import numpy as np
      """
      function will take an input grid, select the values surrounding the input
      coordinates according to the user, then return the array surrounding the coords
      while the values greater than the cutoff distance are set to np.NaN
      this is so plotting is easier
      """
      conversion = self.resolution
      around = around/conversion # must convert A to grids
      rism_grid = np.copy(rism_grid)
      x = np.arange(0,rism_grid.shape[0])
      y = np.arange(0,rism_grid.shape[1])
      z = np.arange(0,rism_grid.shape[2])
      #index_grid = np.ndarray(shape=(rism_grid.shape[0],rism_grid.shape[1],rism_grid.shape[2]))
      index_grid = np.full((rism_grid.shape[0],rism_grid.shape[1],rism_grid.shape[2]),False,dtype=bool)
      #print(index_grid.shape)
      for coords in in_coords:
        #print(coords)
        #x_r,y_r,z_r = int(coords[0]),int(coords[1]),int(coords[2])
        x_r,y_r,z_r = float(coords[0])/conversion,float(coords[1])/conversion,float(coords[2])/conversion
        #mask = (x[np.newaxis,:,:]-x_r)**2+(y[:,np.newaxis,:]-y_r)**2+(z[:,:,np.newaxis]-z_r)**2 < around**2
        #print(x[np.newaxis,:,:])
        index_grid[(x[:,np.newaxis,np.newaxis]-x_r)**2+(y[np.newaxis,:,np.newaxis]-y_r)**2+(z[np.newaxis,np.newaxis,:]-z_r)**2 < around**2] = True
        #rism_grid[(x[np.newaxis,:,:]-x_r)**2+(y[:,np.newaxis,:]-y_r)**2+(z[:,:,np.newaxis]-z_r)**2 < around**2] #= np.nan
        #out_grid = rism_grid[index_grid]#np.where(rism_grid == rism_grid[index_grid])]
      return np.where(index_grid == True,rism_grid,np.nan)

    def extract_grid(self,input_file:str,sele:str='guv',out_name:str='out'):
      """
      This function extracts calc. data from the dump file
      by searching for the value corresponding to the selection string
      This way you dont need to list values to extract them
      """
      assert sele != 'cmd', 'you selected the saved command (cmd)'
      thold = ''.join([chr(i) for i in subprocess.check_output([f"ts4sdump -f {input_file} -l"],shell=True)]).split('\n')
      for val in thold:
        try:
          tmp_str = val.split('@')[0].split()#[1] #.index('guv')
          if tmp_str[1] == sele.strip():
            #print(tmp_str[0])
            value_to_extract = tmp_str[0]
        except IndexError:
          continue
      subprocess.run([f"ts4sdump -f {input_file} -e {value_to_extract} > {out_name}.txt"],shell=True)
          #self.extracted_file = f'{out_name}.txt'
      return

    def select_grid(self,input_string:str='guv',coord_array=None):
      """
      select values in array based on the selection input string
      +++++++++++++++++++++++++++++++++++++++++++
      so far selection string is limited to selection of
      grid-values around single residue names only: i.e. around 4 resname LYS
      But we reccommend using the coordinate array input for more in-depth selections
      :::::: in the future the selection commands will be updated ::::::
      ==================================================
      input_string : value to select and extract
      - default value is guv (atomic density)
      if a selection is made value is passed to -> self.get_coords -> self.extract
      else -> self.extract
      ==================================================
      selection string can include the following key-words
      any output value from 3DRISM, i.e 'guv', 'uuv' or 'coul' etc. etc.
      around 'distance' :: select grid values around coordinates within specified distance (in Angstrom)
      get :: get the exact grid value on the input coordinates (i.e. resname or coord_array)
      resname :: any 3-letter cannonical residue name, upper or lower case
      'laplace' or 'laplacian' or 'lap' or 'grad' or 'gradient' or 'del' :: return laplacian of the whole grid
      this selection can be used in tandem with 'around' and 'get' but will always return the laplacian of the
      original unmasked grid, and NOT the laplacian of the selection region.
      'convolve' or 'log' or 'laplacian of gaussian' :: return the grid convolved with a laplacian of gaussian (LoG)
      filter. Same specification as above
      sigma :: stdev. for gaussian kernel in the LoG filter
      ---------------------------------------------
      ! All other strings will be ignored.
      ! resname will override any array input
      ---------------------------------------------
      example selection
      self.select_grid('LoG of guv with sigma 3 around 4 resname MOL')
      creates laplacian of gaussian of the g(r) grid with a gaussian std. of 3,
      then selects the region around 4 A of all atoms in the residue named 'MOL'
      """
      from types import NoneType
      from numpy import ndarray,array,where,round,append,float64
      lap_flag = False # to laplacian
      conv_flag = False# to convolve
      if type(coord_array) == NoneType:
        coord_flag = False
      elif type(coord_array) == list:
        coord_flag = True
        coord_array = ndarray(coord_array)
        # change array into np.array
      elif type(coord_array) == ndarray:
        coord_flag = True

      # these are the possible names to extract
      # we can add more but for now it is unlikely people
      # will want to select weird values, i.e. hlr
      names = {'guv','uuv','ex','ld','coul'}
      # in the future need to add dictionary so people can
      # use more strings, e.g. select g(r) or local density -> guv, ld
      parser = input_string.lower().split()
      #print(parser)
      item = None # this is our selection value i.e. resname
      for name in names:
        if name in parser:
          val = name
          break # can only select one calculation result

      if ('laplace' or 'laplacian' or 'lap' or 'grad' or 'gradient' or 'del') in parser:
          lap_flag = True
      if ('convolve' or 'log' or 'laplacian of gaussian') in parser:
          conv_flag = True
          if 'sigma' in parser:
              sigma = float(parser[parser.index('sigma')+1])/self.resolution
          else:
              sigma = 1.52/self.resolution # Water VdW radius
      if 'around' in parser:
        #print(parser.index('around'))
        dist = float(parser[parser.index('around')+1])
      else:
        dist = None
      # if no distance is specified then select everything

      if 'get' in parser:
        get_flag = True
        #item = parser[parser.index('resname')+1].upper()
      else:
        get_flag = False
      #print(get_flag)
      if (dist or get_flag) and coord_flag:
        # if the user includes an np.array
        # of coordinates then these will be our coordinates to
        # select around
        item = coord_array

      if 'resname' in parser:
        # will override coordinate array
        # in the future it would be nice to be able to include both
        item = parser[parser.index('resname')+1].upper()
      ############## Now extract ts4s file
      self.extract_grid(f'{self.log}.ts4s',sele=val,out_name=f'{val}_{self.log}')
      # read extracted data into a numpy array
      if lap_flag and (not conv_flag):
          t_grid = self.reader(file_in=f'{val}_{self.log}.txt',laplacian=True)
      if conv_flag and (not lap_flag):
          t_grid = self.reader(file_in=f'{val}_{self.log}.txt',laplacian=False,convolve=True,sigma=sigma)
      else:
          t_grid = self.reader(file_in=f'{val}_{self.log}.txt')#,laplacian=False,file_out='out',dx=False)
      #############

      if (not coord_flag) and (dist):
        # if we have a distance then that means we have a selection
        # so we select values around our coordinates of the selection item
        coords_ = self.select_coords(f'{self.solute_path+self.structure_file}',item)
        #print(coords_)
        return self.select_around(t_grid,coords_,around=dist)

      if (not get_flag) and (not dist):
        # if we dont have a distance that means we dont
        # have a selection and we merely return the entire array
        return t_grid

      if get_flag and (not coord_flag):
        ee = array([],dtype=int)
        #max_x,max_y,max_z = np.float64(t_grid.shape)

        for val in round(self.select_coords(f'{self.solute_path+self.structure_file}',item)/self.resolution):
          ee = append(ee,[int(val[0]),int(val[1]),int(val[2])])

        for max_dim in float64(t_grid.shape):
            if (ee > max_dim).any():
                # if any of our rounded values is the grid limit
                # then round down
                ee[where(ee > max_dim)] = max_dim-1
        if (ee == 0).any():
            # if any of the returned are rounded to zero, then add one grid
            ee[where(ee == 0)] = 1
        # reshape to x,y,z coords
        ee = ee.reshape((int(len(ee)/3),3))
        #return the total array but only where we have selected
        return t_grid[ee[:,0],ee[:,1],ee[:,2]]


      if get_flag and (coord_flag):
        ee = array([],dtype=int)
        #max_x,max_y,max_z = np.float64(t_grid.shape)

        for val in round(item/self.resolution):
          ee = append(ee,[int(val[0]),int(val[1]),int(val[2])])

        for max_dim in float64(t_grid.shape):
            if (ee > max_dim).any():
                # if any of our rounded values is the grid limit
                # then round down
                ee[where(ee > max_dim)] = max_dim-1
        if (ee == 0).any():
            # if any of the returned are rounded to zero, then add one grid
            ee[where(ee == 0)] = 1

        ee = ee.reshape((int(len(ee)/3),3))
        #x[ee]
        return t_grid[ee[:,0],ee[:,1],ee[:,2]]
      else:
        # if we have a distance and coord_array then we
        # return the grid selected around the coord_array
        return self.select_around(t_grid,item,around=dist)
