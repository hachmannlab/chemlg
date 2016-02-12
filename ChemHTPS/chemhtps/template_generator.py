#!/usr/bin/env python

_MODULE_NAME = "template_generator"
_MODULE_VERSION = "v0.1.0"
_REVISION_DATE = "2016-02-04"
_AUTHORS = "Johannes Hachmann (hachmann@buffalo.edu) and William Evangelista (wevangel@buffalo.edu)"
_DESCRIPTION = "This module generates the job template to be used in generating jobs."

# Version history timeline:
# v0.0.1 (2016-01-29): basic implementation
# v0.0.2 (2016-02-02): working with messy nested while loops
# v0.0.3 (2016-02-03): got rid of nested while loops 
# v0.1.0 (2016-02-04): alpha version with full menu functionality

####################################################################################################
# TASKS OF THIS MODULE:
# -genertae the proper input template for a job based on user input
####################################################################################################

####################################################################################################
# TODO:
# possibly add more menus and options
####################################################################################################

import os
import fnmatch
import curses


def runmenu(menu, menu_list):
    """(runmenu):
        This function runs a menu based on a given list.
        :param menu: The curses window
        :param menu_list: The list of options that makes the menu
    """

    curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
    hlight = curses.color_pair(1)
    nohlight = curses.A_NORMAL

    length = len(menu_list)
    pos = 2
    oldpos = None
    keypress = None
    menu.clear()

    while keypress != ord('\n'):
        if pos != oldpos:
            oldpos = pos

            for i in range(length):
                light = nohlight
                if pos == i:
                    light = hlight
                if i == 0:
                    menu.addstr(2, 2, menu_list[0], curses.A_STANDOUT)  # menu title
                elif i == 1:
                    menu.addstr(4, 2, menu_list[1], curses.A_BOLD)  # menu subtitle
                elif 0 != i != 1:
                    menu.addstr(5 + i, 4, "%d - %s" % (i - 1, menu_list[i]), light)

            light = nohlight
            if pos == length - 1:
                light = hlight
            menu.addstr(4 + length, 4, "%d - %s" % (length - 2, menu_list[length - 1]), light)
            menu.refresh()

        keypress = menu.getch()

        if ord('1') <= keypress <= ord(str(length - 2)):
            pos = keypress - ord('0') + 1
        elif keypress == curses.KEY_DOWN:
            if pos < length - 1:
                pos += 1
            else:
                pos = 2
        elif keypress == curses.KEY_UP:
            if pos > 2:
                pos += -1
            else:
                pos = length - 1

    return pos


def showresult(menu, menu_names, options):
    """(showresult):
        This function shows the end results for confirmation
        :param menu: The curses window
        :param menu_names: The names of the options
        :param options: The selected options
    """

    curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
    hlight = curses.color_pair(1)
    nohlight = curses.A_NORMAL

    length = len(options)
    pos = 0
    oldpos = None
    keypress = None
    menu.clear()

    while keypress != ord('\n'):
        if pos != oldpos:
            oldpos = pos
            menu.addstr(2, 2, 'These are the options you have selected', curses.A_STANDOUT)
            menu.addstr(4, 2, 'Select "Yes" to proceed or "No" to restart the process', curses.A_BOLD)

            for i in range(length):
                menu.addstr(6 + i, 4, "%d - %s" % (i + 1, str(menu_names[i]) + ' = ' + options[i]), nohlight)

            yn = ['Yes', 'No']
            for i in range(2):
                light = nohlight
                if pos == i:
                    light = hlight
                menu.addstr(7 + length + i, 4, yn[i], light)
            menu.refresh()
        keypress = menu.getch()

        if keypress == curses.KEY_DOWN:
            if pos < 1:
                pos += 1
            else:
                pos = 0
        elif keypress == curses.KEY_UP:
            if pos > 0:
                pos += -1
            else:
                pos = 1

    return pos


def generate_template():
    """(generate_template):
        This function generates the template to be used for the job_generator.
    """

    cwd = os.getcwd()

    menu = curses.initscr()
    curses.noecho()
    curses.cbreak()
    curses.curs_set(0)
    curses.start_color()
    menu.keypad(1)

    # First select which program will be used eg. ORCA, QChem, etc
    programs = ['Supported Programs', 'Please choose a program', 'ORCA', 'QCHEM', 'Exit Program']
    templates = []
    job_type = ['Supported Job Types', 'Please select a job type', 'Single Point', 'Geometry Optimization',
                'Previous Menu']
    dft_method = ['Supported DFT Methods', 'Please choose a DFT method', 'B3LYP', 'PBE0', 'Previous Menu']
    basis_set = ['Supported Basis Sets', 'Please select a basis set', 'Def2-SVP', 'Def2-TZVP', 'Previous Menu']
    convergence = ['Supported Convergence Levels', 'Please choose a convergence level', 'NormalSCF', 'LooseSCF',
                   'SloppySCF', 'StrongSCF', 'TightSCF', 'VerytightSCF', 'Previous Menu']

    # dictionary of the menus for navigation purposes
    menus = {0: programs, 1: templates, 2: job_type, 3: dft_method, 4: basis_set, 5: convergence}
    menu_names = {0: 'Program', 1: 'Template', 2: 'Job', 3: 'Method', 4: 'Basis', 5: 'Convergence'}
    options = {0: '', 1: '', 2: '', 3: '', 4: '', 5: ''}
    i = 0
    while i <= len(menus):
        if i == len(menus):
            pos = showresult(menu, menu_names, options)
            if pos == 1:
                menus[1] = []
                i = 0
            else:
                break
        pos = runmenu(menu, menus[i])
        options[i] = menus[i][pos]
        if i == 0 and pos == menus[0].index('Exit Program'):
            break
        elif i == 0 and pos != menus[0].index('Exit Program'):
            menus[1].extend(
                ('Available %s Templates' % options[0], 'Please choose a template or select "Make New Template"'))
            for root, directories, filenames in os.walk(cwd + '/job_templates/%s' % options[0]):
                for filename in fnmatch.filter(filenames, '*.inp'):
                    menus[1].append(os.path.join(root, filename))
            menus[1].extend(('Make New Template', 'Previous Menu'))
            i += 1
        elif i != 0 and pos == menus[i].index('Previous Menu'):
            if i == 1:
                menus[1] = []
            i += -1
        elif i == 1 and pos != menus[1].index('Make New Template'):
            i = 6
        else:
            i += 1

    curses.nocbreak()
    curses.echo()
    curses.curs_set(1)
    menu.keypad(0)
    curses.endwin()

    # Now need to either use the selected template or generate a new one
    input_str = ''
    if options[1] != 'Make New Template':  # means an already existing template was selected
        input_str = options[1]
        # print options[1]
    else:  # means a new template was selected
        input_str = cwd
        input_str += '/job_templates/%s' % (options[0] + '/' + options[3] + '_' + options[4] + '.inp')
        with open(input_str, 'w') as input_file:
            tmp = '!DFT ' + options[3] + ' ' + options[4] + ' smallPRINT PRINTBASIS ' + options[5]
            if options[2] == 'Geometry Optimization':
                tmp += ' opt'
            tmp += '\n'
            input_file.write(tmp)
        # for option in options:
            # print options[option]

    return input_str
