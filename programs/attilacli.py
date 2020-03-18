#!/usr/bin/env python3

"""

attila_cli.py - Command line interface to setup Attila projects

Usage:
    python3 attila_cli.py # TODO: finish command line parser

Copyright: (c) 2019-2020 Matheus Cardoso <github.com/cardosaum>

License: Apache 2.0 <apache.org/licenses/LICENSE-2.0>

"""


import argparse
import subprocess
import pathlib
import os
import readline
import pprint
import collections
import time

#############################################
# setup variables to be used on the hole program
readline.set_completer_delims(" \t\n=")
readline.parse_and_bind("tab: complete")
yes_list = ['y']
no_list = ['n', '']
str_only_yes_or_no = 'Please, enter "Y" for "Yes" or "n" for "No".'

# documentation string
# TODO: does not have a better way to do this?
# maybe using `argparse`, `click` or `optparse` libraries
attila_commands = {
    'commands:':
    {
        'CTRL-C': 'quit ATTILA; abort analysis',
        'TAB': 'autocomplete a path'
    },
    'configuration parameters:':
    {
        'Configuration files exist (y or n)': 'type \'y\' if you already have configuration files\n\
                    type \'n\' or press ENTER key if you prefer to let ATTILA create the conf{}iguration files'.format((' '*64)),
        'Path of the configuration file of VH libraries': 'location of the configuration file of VH libraries',
        'Path of the configuration file of VL libraries': 'location of the configuration file of VL libraries',
        'Project Name': 'name of the directory that will be created by ATTILA to save output {}files'.format((' '*67)),
        'Directory to save project': 'the directory where the project will be saved',
        'Reads are paired-end (y or n)': 'type \'y\' or press ENTER key for yes; type \'n\' for no',
        'Minimum read length': 'default value is 300 pb; type \'y\' to change default; type \'n\' or press {}ENTER key to use default value\n\
                    if you choose to change default value, the new read length must be an i{}nteger number'.format((' '*64), (' '*64)),
        'Minimum base quality': 'default value is 20; type \'y\' to change default; type \'n\' or press ENTE{}R key to use default value\n\
                    if you choose to change default value, the new base quality must be an {}integer number'.format((' '*64), (' '*64)),
        'Number of candidates to rank': 'number of candidate clones that ATTILA will try to find in VH and VL li{}braries\n\
                    the number must be an integer'.format((' '*64))
    },
    'Parameters for paired-end reads:':
    {
        'Path of fastq file of VH R0 reads r1': 'location of the fastq file containing reads r1 from initial VH library',
        'Path of fastq file of VH R0 reads r2': 'location of the fastq file containing reads r2 from initial VH library',
        'Path of fastq file of VH RN reads r1': 'location of the fastq file containing reads r1 from final VH library',
        'Path of fastq file of VH RN reads r2': 'location of the fastq file containing reads r2 from final VH library',
        'Path of fastq file of VL R0 reads r1': 'location of the fastq file containing reads r1 from initial VL library',
        'Path of fastq file of VL R0 reads r2': 'location of the fastq file containing reads r2 from initial VL library',
        'Path of fastq file of VL RN reads r1': 'location of the fastq file containing reads r1 from final VL library',
        'Path of fastq file of VL RN reads r2': 'location of the fastq file containing reads r2 from final VL library'
    },
    'Parameters for single-end reads:':
    {
        'Path of fastq file of VH R0': 'location of fastq file containing reads from initial VH library',
        'Path of fastq file of VH RN': 'location of fastq file containing reads from initial VH library',
        'Path of fastq file of VL R0': 'location of fastq file containing reads from initial VH library',
        'Path of fastq file of VL RN': 'location of fastq file containing reads from initial VH library'
    }
}

# `yes or no` related
vals = collections.defaultdict(dict)
vals['yes_list'] = ['y']
vals['no_list'] = ['n', ''] # empty string will be considered as "No"
vals['yn_valid_inputs'] = vals['yes_list'] + vals['no_list']
def hint():
    """stores a hint of type '[y/N]' (capital letter depending if empty string is in `yes_list` or `no_list`)"""

    hint = '['
    if '' in vals['yes_list']:
        hint += 'Y/n]'
    else:
        hint += 'y/N]'

    return hint
vals['yn_hint'] = hint()

# string related
# set default character separator
vals['default_char_sep'] = '*'
vals['default_prompt'] = '> '

def config_default_file_names():
    """set default configuration file names"""
    for vX in ['VH', 'VL']:
        vals[f'default_config_file_{vX.lower()}'] = pathlib.Path(f"{utils_config_get_settings_value(1)}_{vX}.cfg")

def config_default_attila_program_names():
    """set default attila programs' names"""
    vals['attila_programs'] = ['autoiganalysis3.pl', 'translateab9', 'frequency_counter3.pl', 'find_duplicates7.pl', 'get_nsequences.pl', 'numberab.pl', 'convertofasta.pl', 'get_ntsequence2.pl', 'rscript_creator.pl', 'html_creator.pl', 'parserid.pl', 'statscript_creator.pl']

config_default_attila_program_names()
# attila settings
    # settings[1]      Project name
    # settings[2]      Project path
    # settings[3]      Attila package path
    # settings[4]      Reads are paired-end (0/1)
    # settings[5]      VH R0 reads r1 path
    # settings[6]      VH R0 reads r2 path
    # settings[7]      VH RN reads r1 path
    # settings[8]      VH RN reads r2 path
    # settings[9]      VL R0 reads r1 path
    # settings[10]      VL R0 reads r2 path
    # settings[11]      VL RN reads r1 path
    # settings[12]      VL RN reads r2 path
    # settings[13]      VH R0 path
    # settings[14]      VH RN path
    # settings[15]      VL R0 path
    # settings[16]      VL RN path
    # settings[17]      IgBlast package path
    # settings[18]      Minimum read length
    # settings[19]      Minimum base quality
    # settings[20]      Number of candidates to rank
    # settings[21]      VH file configuration path
    # settings[22]      VL file configuration path
s = collections.defaultdict(dict)
def populate_settings(num=22):
    """Populate settings variables until setting number `num`
Output expected:

s = {
    '1': ['Project name', ''],
    '2': ['Project path', ''],
    '3': ['Attila package path', ''],
    '4': ['Reads are paired-end (True/False)', ''],
    '5': ['VH R0 reads r1 path', ''],
    '6': ['VH R0 reads r2 path', ''],
    '7': ['VH RN reads r1 path', ''],
    '8': ['VH RN reads r2 path', ''],
    '9': ['VL R0 reads r1 path', ''],
    '10': ['VL R0 reads r2 path', ''],
    '11': ['VL RN reads r1 path', ''],
    '12': ['VL RN reads r2 path', ''],
    '13': ['VH R0 path', ''],
    '14': ['VH RN path', ''],
    '15': ['VL R0 path', ''],
    '16': ['VL RN path', ''],
    '17': ['IgBlast package path', ''],
    '18': ['Minimum read length', ''],
    '19': ['Minimum base quality', ''],
    '20': ['Number of candidates to rank', ''],
    '21': ['VH file configuration path', ''],
    '22': ['VL file configuration path', '']
    }
    """
    settings_names = ['Project name', 'Project path', 'Attila package path', 'Reads are paired-end (True/False)', 'VH R0 reads r1 path', 'VH R0 reads r2 path', 'VH RN reads r1 path', 'VH RN reads r2 path', 'VL R0 reads r1 path', 'VL R0 reads r2 path', 'VL RN reads r1 path', 'VL RN reads r2 path', 'VH R0 path', 'VH RN path', 'VL R0 path', 'VL RN path', 'IgBlast package path', 'Minimum read length', 'Minimum base quality', 'Number of candidates to rank', 'VH file configuration path', 'VL file configuration path']
    for i in range(1, num+1):
        s.setdefault(str(i), [settings_names[i-1], ''])
    return s
populate_settings()
#############################################



# TODO: simple function to interact with user


# TODO: see if configuration already exists

# TODO: set all configuration by asking user
    # TODO: handle all paths correctly, expanding them

# TODO: show user all config and ask if correct

# TODO: once confirmed that config is right, create all
# needed directories and config files





def update_vals():
    """ Get values for terminal """

    # terminal related
    info = os.get_terminal_size()
    vals['terminal_info'] = info
    vals['terminal_size'] = info.columns


    return vals


def user_greetings(text="ATTILA: Automated Tool for Immunoglobulin Analysis", clear_screen=True):
    """Greeting user"""
    
    # to display formated text correctly
    update_vals()

    if clear_screen:
        subprocess.Popen('clear')
        time.sleep(0.01) # some terminals does not show above text if `clear` and `print` execute
                        # at near the same time

    # Greeting user
    print(vals['default_char_sep'] * vals['terminal_size'])
    print(text.center(vals['terminal_size']))
    print(vals['default_char_sep'] * vals['terminal_size'])
    print()

def config_check_existing():
    """Search recursively current working directory for existing configuration files"""
    pass
    # TODO: write this function. 
    #       *AND* ensure this does not slows down the program
    #       (this function shoul have nearly instant execution time!)


def utils_str_convert_yn_tf(yn:str):
    """Convert 'y' or 'n' to `True` or `False`"""
    yn = yn.lower().strip()
    # verify if `yn` is an acceptable input
    err = True
    for i in vals['yn_valid_inputs']:
        i = i.lower()
        if i == yn:
            err = False
    assert not err, ValueError('Input must be in `yn_valid_inputs`')

    if yn in vals['yes_list']:
        return True
    elif yn in vals['no_list']:
        return False
    else:
        raise ValueError(f'Value {yn} was not expected')

def utils_str_path_absolute(path:str):
    """Return absolute path"""
    if path.startswith('~'):
        path = path.replace('~', str(pathlib.Path.home()), 1)
    p = pathlib.Path(path).resolve()
    return p

def utils_str_length_bigger(st:collections.defaultdict):
    """Return the longest string length in a given list"""
    m = 0
    for i in st.items():
        i = i[1][0]
        if len(i) > m:
            m = len(i)
    return m

def utils_settings_get_number_of_settings(st:collections.defaultdict):
    """Return the maximun value of specified settings in `settings` dictionary"""
    m = max([ int(i[0]) for i in st.items() ])
    return m

def user_ask_yn(text):
    """Ask user a simple 'Yes' or 'No' question"""

    while True:
        print(f"{text} {vals['yn_hint']}")
        answer = input(vals['default_prompt'])
        if answer in vals['yn_valid_inputs']:
            break
        
        print('Please, you must insert only:\n"Y" for "Yes";\n"N" for "No".\n')


    return utils_str_convert_yn_tf(answer)

def user_ask_config():
    """Ask user if there already exist configuration files"""
    exist_config = user_ask_yn('Configuration files exist?')
    return exist_config

def user_ask_setting(text:str, num:int, path=False, empty=False, file=False, fastq=False, default=None, inputInt=False):
    """Ask user for setting number `num`"""

    # path can't be a empty string `""`
    if path:
        empty = False

    # file must have a validated path
    if file:
        path = True

    # set valid fastq files
    if fastq:
        fastq = ['.fq', '.fastq']

    # if `default` is true, empty inputs are allowed.
    # also, format `text` to show this `default` value
    if default:
        empty = True
        text = text.strip() + f"\n(default={default})"

    while True:
        print(text)
        answer = input(vals['default_prompt'])

        # if empty is `False`, tell it for user and ask again
        if not empty and not answer:
            print('This configuration can *not* be empty. Please, insert a valid value.\n')
            continue
        if inputInt:
            if default and not answer:
                answer = default
            try:
                answer = int(answer)
                if answer <= 0:
                    raise ValueError('This value must be a positive number')
            except:
                print('This configuration must be a *Positive Integer Number*.\n')
                continue
        if default and not answer:
            answer = default
        if path:
            path = pathlib.Path(answer).expanduser().resolve()
            if not file and not path.is_dir():
                print(f"Path {path} is not a regular directory.\nPlease, choose a valid one.\n\n")
                mkdir_yn = user_ask_yn(f'Do you want to create path:\n· "{path}"?')
                if mkdir_yn:
                    path.mkdir(parents=True)
                    s[str(num)][-1] = answer
                    break
                continue
            elif file and not path.is_file():
                print(f"File {path} is not a regular file.\nPlease, choose a valid one.\n\n")
                continue
            elif file and path.suffix not in fastq:
                print(f"File {path} is a regular file, but it's extension is *not* a valid one.\nPlease, choose a file that have one of this valid extensions: \"{', '.join(fastq)}\"")
            else:
                answer = str(path)

        s[str(num)][-1] = answer
        break
    return s[str(num)]

def user_ask_reads_paired_end(paired:bool):
    """Aks user for fastq files for all reads"""

    # counter for `setting` number
    if paired:
        num = 5
    else:
        num = 13

    for vX in ['H', 'L']:
        for rX in ['0', 'N']:
            if paired:
                for readX in ['1', '2']:
                    user_ask_setting(f"Enter path for fastq file of V{vX} R{rX} read r{readX}", num, file=True, fastq=True)
                    num += 1
            else:
                user_ask_setting(f"Enter path for fastq file of V{vX} R{rX}", num, file=True, fastq=True)
                num += 1


def user_show_configs_all(s:collections.defaultdict, fill_character='_', text=' Current User Configuration '):
    """Print all current configs to user"""

    # For `settings` number and value, print it's value to user
    l = utils_str_length_bigger(s)
    m = 0
    final_text = ''
    # This first loop is for get max length needed
    for n, v in s.items():
        t1 = f"· ({n}) {v[0]}:\t ".ljust(l)
        # t2 = f"\t {v[-1]}"
        # t3 = t1 + t2
        if len(t1) > m:
            m = len(t1)

    print()
    print(vals['default_char_sep'] * vals['terminal_size'])
    print(text.center(vals['terminal_size']))
    print(vals['default_char_sep'] * vals['terminal_size'])
    print()
    # This loop actualy print the correct message to user
    for n, v in s.items():
        t1 = f"· ({n}) {v[0]}: ".ljust(m, fill_character)
        t2 = f" {v[-1]}"
        t3 = t1 + t2
        print(t3)
    print()

def user_ask_config_change_one():
    """Ask user to change one specific configuration"""

    g = [int(i) for i in s.keys()]
    while True:
        print()
        print('Which configuration do you want to change?')
        try:
            c_num = int(input('> ').strip())
        except:
            print(f'Please, you must insert a *Valid Integer Number*.\nThe valid ones are in range "{min(g)}-{max(g)}"\n')
            continue
        if c_num < min(g) or c_num > max(g):
            print(f'Please, you must insert a *Valid Integer Number*.\nThe valid ones are in range "{min(g)}-{max(g)}"\n')

        else:
            break
    return c_num


def user_ask_configs_all():
    """Show user all current configs and ask if all of them are correct.
    If they are not correct, show a menu to let user change the wrong ones.
    """
    global s

    while True:
        user_show_configs_all(s)
        change = user_ask_yn('Do you want to change any configuration?')
        if not change:
            break
        setting_num = user_ask_config_change_one()
        user_ask_default_settings(setting_num)


def utils_config_get_attila_packge_path():
    """Get path where attila is saved"""

    p = list(pathlib.Path(__file__).parents)
    for i in p:
        if str(i).endswith(str(pathlib.Path.joinpath(pathlib.Path('attila'), 'programs'))):
            i = i.parent
            return i

    # If we couldn't find the path, raise a error
    raise FileNotFoundError('Attila package directory could not be found.')

def user_ask_default_settings(n, manualy=True):
    """Set default settings to all user questions"""
    if n == 1:
        user_ask_setting("Enter project name:", 1, empty=False)
    elif n == 2:
        user_ask_setting("Enter directory to save the project:", 2, path=True)
    elif n == 3:
        s['3'][-1] = utils_config_get_attila_packge_path()
    elif n == 4:
        s['4'][-1] = user_ask_yn("Reads are paired-end?") 
        user_ask_reads_paired_end(s['4'][-1])
    elif n == 18:
        user_ask_setting("Minimum read length:", 18, default=300)
    elif n == 19:
        user_ask_setting("Minimum base quality:", 19, default=20)
    elif n == 20:
        user_ask_setting("Enter number of candidates to rank:", 20, inputInt=True)
    else:
        if manualy:
            print(f'Sorry, but the setting number {n} is not intended to be configured manualy.')

def utils_config_get_settings_value(n:int):
    return s[str(n)][-1]


def config_create_file(file:pathlib.Path):
    """Create configuration files"""

    # TODO: This file can be way better if we use a more standard format of file
    #       such as JSON. I only wrote this function for backward compatibility
    #       reasons.
    #       Can we use a more suitable format in later versions?
    #       - Matheus Cardoso , 12/03/2020.

    # set a helper variable
    config_default_file_names()

    # just an alias
    gv = utils_config_get_settings_value

    # create file for VH and VL reads
    for f in ['VH', 'VL']:
        if f == 'VH':
            t = 0
        else:
            t = 1
        if gv(4):
            gv4 = 1
        else:
            gv4 = 0
        text = []
        text.append('-' * 70)
        text.append("# [ Section: files and directories ]")
        text.append(f"projectname: {gv(1)}")
        text.append(f"projectdir: {gv(2)}")
        text.append(f"packagedir: {gv(3)}")
        text.append(f"igblastdir: {gv(17)}")
        if f == 'VH' and gv(4):
            text.append(f"input1r1dir: {gv(5)}")
            text.append(f"input1r2dir: {gv(6)}")
            text.append(f"input2r1dir: {gv(7)}")
            text.append(f"input2r2dir: {gv(8)}")
        if f == 'VH' and not gv(4):
            text.append(f"input1dir: {gv(13)}")
            text.append(f"input2dir: {gv(14)}")
        if f == 'VL' and gv(4):
            text.append(f"input1r1dir: {gv(9)}")
            text.append(f"input1r2dir: {gv(10)}")
            text.append(f"input2r1dir: {gv(11)}")
            text.append(f"input2r2dir: {gv(12)}")
        if f == 'VL' and not gv(4):
            text.append(f"input1dir: {gv(15)}")
            text.append(f"input2dir: {gv(16)}")
        text.append("# [ Section: analysis arguments ]")
        text.append(f"libtype: {t}")
        text.append(f"listsize: {gv(20)}")
        text.append(f"pairedend: {gv4}")
        text.append(f"minlen: {gv(18)}")
        text.append(f"minqual: {gv(19)}")
        final_text = '\n'.join(text)
        with vals[f'default_config_file_{f.lower()}'].open('w') as c:
            c.write(final_text)


def config_create_symlinks(path:pathlib.Path, dir_symlink="ATTILASymLinks"):
    """create all needed symlinks for project.
    `path` correspond to `s['2'][-1]`. i.e. `Project path` specified by the user
    """

    # first, let's figure out if the Symlink directory already exists
    if path.exists():
        # TODO: what we are supposed to do if this directory already exists?
        symlink_path = pathlib.Path.joinpath(path, dir_symlink)
        pathlib.Path.mkdir(symlink_path, parents=True, exist_ok=True)

        for p in vals['attila_programs']:
            source_path = pathlib.Path.joinpath(utils_config_get_settings_value(3), 'programs', p)
            dest_path = pathlib.Path.joinpath(symlink_path, p)
            try:
                dest_path.symlink_to(source_path)
            except FileExistsError:
                pass

        # raise FileExistsError('ATTILA is not yet configured to handle a Symlink Directory that already exists.\n\nPlease contact either:\n\n· Waldeyr (https://github.com/waldeyr)\n\nor\n\n· Matheus Cardoso (https://github.com/cardosaum)')

    # If symlinks' directory does not exists, we'll create it, with all needed files
    else:

        symlink_path = pathlib.Path.joinpath(path, dir_symlink)
        pathlib.Path.mkdir(symlink_path, parents=True)

        for p in vals['attila_programs']:
            source_path = pathlib.Path.joinpath(utils_config_get_settings_value(3), 'programs', p)
            dest_path = pathlib.Path.joinpath(symlink_path, p)
            dest_path.symlink_to(source_path)


def config_run_analisys():
    """create needed directories and run analysis"""

    # subprocess.Popen('clear', shell=True)

    # First, let's create the needed folders
    print('Creating project directory')

    jp = pathlib.Path.joinpath
    directories = []
    project_path = jp(pathlib.Path(utils_config_get_settings_value(2)), utils_config_get_settings_value(1))
    report_path = jp(project_path, 'Report')

    directories.append(project_path)
    directories.append(report_path)

    [ pathlib.Path.mkdir(i, parents=True) for i in directories ]

    # Now, we run the analisys
    for vX in ['VH', 'VL']:

        print(f'Running {vX} analysis ...')
    
        vx_error_log = pathlib.Path.joinpath(project_path, f'{vX.lower()}error.log')
        symlink_path = pathlib.Path.joinpath(pathlib.Path(utils_config_get_settings_value(2)), 'ATTILASymLinks')
        perl_file = pathlib.Path.joinpath(symlink_path, 'autoiganalysis3.pl')
        subprocess.Popen(f'time perl {perl_file} {vx_error_log}', shell=True)

        print('-' * vals['terminal_size'])
        print(f'{vX} Analysis Completed')
        print('-' * vals['terminal_size'])




def utils_files_find(path:pathlib.Path, pattern, num=1):
    """find files in `path` by regex rules specified by `pattern`""" 

    files = sorted(path.rglob(pattern))
    return files

def ff(path:pathlib.Path, pattern):
    """An alias function for `utils_files_find`"""
    return utils_files_find(path, pattern)


def utils_find_vx(files:list, vx:str):
    """find file belonging only to vX, where vX can assume either VH or VL"""
    vx = vx.upper()
    for i in files:
        f = str(i).upper()
        if vx in f:
            pass
            # print(i)
        else:
            print(i)




def config_create_web_page():
    """create the web page containing the results"""

    project_path = pathlib.Path.joinpath(pathlib.Path(utils_config_get_settings_value(2)), utils_config_get_settings_value(1))
    report_path = pathlib.Path.joinpath(project_path, 'Report')

    project_path.mkdir(exist_ok=True)
    report_path.mkdir(exist_ok=True, parents=True)

    # a dictionary to store all nedded files
    f = collections.defaultdict(dict)

    for vX in ["VH", "VL"]:
        f[f'numbered_{vX}'] = ff(project_path, f"*{vX}/*numbered.fasta")
        f[f'germline_{vX}'] = ff(project_path, f"*{vX}/*germlineclassification.txt")
        f[f'plot1_{vX}'] = ff(project_path, f"*length_{vX.lower()}.png")
        f[f'plot2_{vX}'] = ff(project_path, f"*task_{vX.lower()}.png")

    command = []
    command.append("perl")
    command.append(str(pathlib.Path.joinpath(pathlib.Path(utils_config_get_settings_value(2)), 'ATTILASymLinks', 'html_creator.pl')))
    command = ' '.join(command)

    command = f"{command} {f['numbered_VH']} {f['numbered_VL']} {f['germline_VH']} {f['germline_VL']} {report_path}/Report.html {plot1_VH} {plot2_VH} {plot1_VL} {plot2_VL} {pathlib.Path.joinpath(project, 'VH', 'vhSequenceCounting.csv')} {pathlib.Path.joinpath(project, 'VL', 'vlSequenceCounting.csv')} {pathlib.Path.joinpath(report_path, 'vhoutputRstats.txt')} {pathlib.Path.joinpath(report_path, 'vloutputRstats.txt')} > {pathlib.Path.joinpath(report_path, 'webpage.log')} 2>&1"
    subprocess.Popen(command)

def config_check_analisys_result():
    """check if analysis suceeded"""

    project_path = pathlib.Path.joinpath(utils_config_get_settings_value(2), utils_config_get_settings_value(1))
    report_path = pathlib.Path.joinpath(project_path, 'Report')
    webpage_file = pathlib.Path.joinpath(report_path, 'webpage.log')
    with webpage_file.open() as f:
        text = f.read()
        if text:
            print('Could not create analysis report')
            return False
        else:
            print('Analysis report is ready!')
            return True



def flow():
    """Flow control for program"""
    user_greetings()

    # Ask user if configutarion files already exist
    if user_ask_config():
        # if configuration does already exists, use it
        user_ask_setting(f"Enter path for configuration file of V{vX} libraries", num)
        pass

    else:
        # TODO: if configuration does not already exists, create it
        for f in range(utils_settings_get_number_of_settings(s)):
            user_ask_default_settings(f, manualy=False)

        # ask user if all inserted configs are indeed correct
        user_ask_configs_all()

        # now, we need to write the configs to a file
        config_create_file(pathlib.Path(utils_config_get_settings_value(1)))

        # and, finaly, we will create the relevant symlinks
        config_create_symlinks(pathlib.Path(utils_config_get_settings_value(2)))

        # Now we only need to run the analisys, create the webpage and check their result
        config_run_analisys()
        config_create_web_page()
        config_check_analisys_result()


    # Now that we have all the configurations, ask user if this all is correct
    # TODO: 
    pass




def main():
    """Execute main program"""
    flow()

def test():
    # user_ask_setting("Enter number of candidates to rank:", 20, inputInt=True)
    # pprint.pprint(s)
    # print(pathlib.Path(__file__).resolve())
    # print(utils_settings_get_number_of_settings(s))
    # pprint.pprint(s)
    # pprint.pprint(s[str(1)][-1])
    # s[str(1)][-1] = 'FOO'
    # pprint.pprint(s)
    # pprint.pprint(s[str(1)][-1])

    # g = ff(pathlib.Path.joinpath(utils_config_get_attila_packge_path()), "*V[HL]*")
    # g = ff(pathlib.Path("/home/matheus/trash_mcs/tm/t"), "*vh/*fasta")
    # pprint.pprint(g)
    # utils_find_vx(g, "vx")

    # s['2'][-1] = pathlib.Path(__file__).parent.parent.parent.parent.parent.parent
    # config_create_web_page()


    # config_create_file(pathlib.Path('arquvs'))
    # user_ask_configs_all()

    pass


if __name__ == '__main__':
    main()
    # test()
