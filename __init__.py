from shutil import which 
names = ['eprism3d','ts4sdump','gmxtop2solute','generate-idc.sh']

################################################################################
kernel_warn = """
it appears you dont have the Episol kernel sourced or installed !
Please see the EPISOL_kernel_install_instrcutions.md for
a basic overview on how to install the kernel and make sure
you have the $PATH set to the kernel (eprism3d)
'Visit https://github.com/EPISOLrelease/EPIPY for more
"""
ts4s_warn = """
We cannot find ts4sdump. it appears you dont have the executable sourced properly
"""
top2sol_warn ="""
We cannot find gmxtop2solute it appears you dont have the executable sourced properly
""" 
gen_idc_warn = """
We cannot find generate-idc.sh it appears you dont have the executable sourced properly
"""
################################################################################
kernel_res = """
You will be unable to run ANY 3DRSIM calculations
"""
ts4s_res = """
You will be unable to extract and read your calculation output files
"""
top2sol_res ="""
You will be unable to convert .top files to .solute files and hence be unable to generate IDC enabled files
"""
gen_idc_res = """
You will be unable to generate ion-dipole correction (IDC) enabled .solute files
"""
################################################################################

warns = [kernel_warn,ts4s_warn,top2sol_warn,gen_idc_warn]
results = [kernel_res,ts4s_res,top2sol_res,gen_idc_res]

# this is kind of a goofy way of doing this but it works

err_dict = dict(zip(names,[which(name) for name in names]))
warn_dict = dict(zip(names,[warn for warn in warns]))
result_dict = dict(zip(names,[res for res in results]))

outcomes = []
for name in names:
    if err_dict[name] is None:
        outcomes.append(result_dict[name])
        print('!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!')
        print(warn_dict[name])
outcomes = []
if outcomes  == []:
    print("All Checks complete\nGood to go!")
else:
    print('!!!!!!!!!!!!!!!!!!! OUTCOMES !!!!!!!!!!!!!!!!!!!')
    for val in outcomes:
        print(val) 

from episol.epipy import epipy
