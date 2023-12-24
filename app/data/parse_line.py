"""
SPCAT .cat file parsing script adapted from molsim package
developed by Brett McGuire et al.
"""
import numpy as np
import re
import scipy.constants
from data.class_parse_catfile import (
    Catalog, Level, PartitionFunction, _read_txt)

# Boltzmann's constants in cm-1/K
kcm = scipy.constants.k * (scipy.constants.h**-1) * \
    ((scipy.constants.c * 100)**-1)
ccm = scipy.constants.c * 100  # speed of light in cm/s


def _make_qnstr(qn1, qn2, qn3, qn4, qn5, qn6, qn7, qn8):
    qn_list = [qn1, qn2, qn3, qn4, qn5, qn6, qn7, qn8]
    tmp_list = [str(x).zfill(2) for x in qn_list if x is not None]
    return ''.join(tmp_list)


def _read_spcat(filein):
    '''
    Reads an SPCAT catalog and returns spliced out numpy arrays.
    Catalog energy units for elow are converted to K from cm-1.
    '''

    # read in the catalog
    raw_arr = _read_txt(filein)

    # set up some basic lists to populate
    frequency, freq_err, logint, dof, elow, gup, tag, qnformat, \
        qn1, qn2, qn3, qn4, qn5, qn6, qn7, qn8, qn9, qn10, qn11, \
        qn12 = ([] for _ in range(20))

    # split everything out
    for x in raw_arr:
        """
        if there is a * in there, then the catalog was
        simulated too high for SPCAT format and we need to skip the line.
        """
        if '*' in x:
            continue
        else:
            frequency.append(x[:13].strip())
            freq_err.append(x[13:21].strip())
            logint.append(x[21:29].strip())
            dof.append(x[29:31].strip())
            elow.append(x[31:41].strip())
            gup.append(x[41:44].strip())
            tag.append(x[44:51].strip())
            qnformat.append(x[51:55].strip())
            qn1.append(x[55:57].strip())
            qn2.append(x[57:59].strip())
            qn3.append(x[59:61].strip())
            qn4.append(x[61:63].strip())
            qn5.append(x[63:65].strip())
            qn6.append(x[65:67].strip())
            qn7.append(x[67:69].strip())
            qn8.append(x[69:71].strip())
            qn9.append(x[71:73].strip())
            qn10.append(x[73:75].strip())
            qn11.append(x[75:77].strip())
            qn12.append(x[77:].strip())

    # now go through and fix everything into the appropriate
    # formats and make numpy arrays as needed

    # we start with the easy ones that don't have nonsense letters
    frequency = np.array(frequency)
    frequency = frequency.astype(np.float64)
    freq_err = np.array(freq_err)
    freq_err = freq_err.astype(np.float64)
    logint = np.array(logint)
    logint = logint.astype(np.float64)
    dof = np.array(dof)
    dof = dof.astype(np.int64)
    elow = np.array(elow)
    elow = elow.astype(np.float64)
    tag = np.array(tag)
    tag = tag.astype(np.int64)
    qnformat = np.array(qnformat)
    qnformat = qnformat.astype(np.int64)

    # convert elow to Kelvin

    elow /= kcm

    # now we use a sub-function to fix the letters and
    # +/- that show up in gup, and the qns, and returns floats
    # Edited by Jasmine to account for upper and lower case alphabets
    # that represent positive and negative quantum numbers respectively.
    def _fix_spcat(x):
        '''Fixes letters in something that's read in and returns floats'''

        # fix blanks - we just want them to be nice nones
        # rather than empty strings

        if x == '':
            return None

        sub_dict = {'A': 100,
                    'B': 110,
                    'C': 120,
                    'D': 130,
                    'E': 140,
                    'F': 150,
                    'G': 160,
                    'H': 170,
                    'I': 180,
                    'J': 190,
                    'K': 200,
                    'L': 210,
                    'M': 220,
                    'N': 230,
                    'O': 240,
                    'P': 250,
                    'Q': 260,
                    'R': 270,
                    'S': 280,
                    'T': 290,
                    'U': 300,
                    'V': 310,
                    'W': 320,
                    'X': 330,
                    'Y': 340,
                    'Z': 350,
                    'a': -10,
                    'b': -20,
                    'c': -30,
                    'd': -40,
                    'e': -50,
                    'f': -60,
                    'g': -70,
                    'h': -80,
                    'i': -90,
                    'j': -100,
                    'k': -110,
                    'l': -120,
                    'm': -130,
                    'n': -140,
                    'o': -150,
                    'p': -160,
                    'q': -170,
                    'r': -180,
                    's': -190,
                    't': -200,
                    'u': -210,
                    'v': -220,
                    'w': -230,
                    'x': -240,
                    'y': -250,
                    'z': -260,
                    }
        # spinv said the minimum quantum number is -259, but I
        # think it's a mistake, it should be able to store down to -269.

        # replace any characters that are NOT alphabets with an empty string
        alpha = re.sub('[^a-zA-Z]+', '', x)

        if alpha == '':
            return int(x)
        else:
            alphabet = sub_dict.get(alpha, 0)
            if alphabet < 0:
                # if it's negative, we need to subtract
                return alphabet - int(x[1])
            else:
                # if it's positive, we need to add
                return alphabet + int(x[1])

    # run the other arrays through the fixer, then convert
    # them to what they need to be

    gup = [_fix_spcat(x) for x in gup]
    gup = np.array(gup)
    gup = gup.astype(int)
    qn1 = [_fix_spcat(x) for x in qn1]
    qn1 = np.array(qn1)
    qn1 = qn1.astype(int) if all(y is not None for y in qn1) is True else qn1
    qn2 = [_fix_spcat(x) for x in qn2]
    qn2 = np.array(qn2)
    qn2 = qn2.astype(int) if all(x is not None for x in qn2) is True else qn2
    qn3 = [_fix_spcat(x) for x in qn3]
    qn3 = np.array(qn3)
    qn3 = qn3.astype(int) if all(x is not None for x in qn3) is True else qn3
    qn4 = [_fix_spcat(x) for x in qn4]
    qn4 = np.array(qn4)
    qn4 = qn4.astype(int) if all(x is not None for x in qn4) is True else qn4
    qn5 = [_fix_spcat(x) for x in qn5]
    qn5 = np.array(qn5)
    qn5 = qn5.astype(int) if all(x is not None for x in qn5) is True else qn5
    qn6 = [_fix_spcat(x) for x in qn6]
    qn6 = np.array(qn6)
    qn6 = qn6.astype(int) if all(x is not None for x in qn6) is True else qn6
    qn7 = [_fix_spcat(x) for x in qn7]
    qn7 = np.array(qn7)
    qn7 = qn7.astype(int) if all(x is not None for x in qn7) is True else qn7
    qn8 = [_fix_spcat(x) for x in qn8]
    qn8 = np.array(qn8)
    qn8 = qn8.astype(int) if all(x is not None for x in qn8) is True else qn8
    qn9 = [_fix_spcat(x) for x in qn9]
    qn9 = np.array(qn9)
    qn9 = qn9.astype(int) if all(x is not None for x in qn9) is True else qn9
    qn10 = [_fix_spcat(x) for x in qn10]
    qn10 = np.array(qn10)
    qn10 = qn10.astype(int) if all(
        x is not None for x in qn10) is True else qn10
    qn11 = [_fix_spcat(x) for x in qn11]
    qn11 = np.array(qn11)
    qn11 = qn11.astype(int) if all(
        x is not None for x in qn11) is True else qn11
    qn12 = [_fix_spcat(x) for x in qn12]
    qn12 = np.array(qn12)
    qn12 = qn12.astype(int) if all(
        x is not None for x in qn12) is True else qn12

    # make the qnstrings
    qn_list_up = [_make_qnstr(qn1, qn2, qn3, qn4, qn5, qn6, None, None)
                  for qn1, qn2, qn3, qn4, qn5, qn6 in zip(
                      qn1, qn2, qn3, qn4, qn5, qn6)]
    qn_list_low = [_make_qnstr(qn7, qn8, qn9, qn10, qn11, qn12, None, None)
                   for qn7, qn8, qn9, qn10, qn11, qn12 in zip(
                       qn7, qn8, qn9, qn10, qn11, qn12)]

    split_cat = {
        'frequency': 	frequency,
        'freq_err':	freq_err,
        'logint':	logint,
        'dof':	dof,
        'elow':	elow,
        'gup':	gup,
        'tag':	tag,
        'qnformat':	qnformat,
        'qn1up':	qn1,
        'qn2up':	qn2,
        'qn3up':	qn3,
        'qn4up':	qn4,
        'qn5up':	qn5,
        'qn6up':	qn6,
        'qn7up':	np.full(len(frequency), None),
        'qn8up':	np.full(len(frequency), None),
        'qnup_str':	np.array(qn_list_up),
        'qn1low':	qn7,
        'qn2low':	qn8,
        'qn3low':	qn9,
        'qn4low':	qn10,
        'qn5low':	qn11,
        'qn6low':	qn12,
        'qn7low':	np.full(len(frequency), None),
        'qn8low':	np.full(len(frequency), None),
        'qnlow_str':	np.array(qn_list_low),
        'notes':	f'Loaded from file {filein}'
    }

    return split_cat


def _load_catalog(filein):
    '''
    Reads in a catalog file of the specified type and returns a catalog
    object. Optionally accepts a catdict dictionary to preload the
    catalog object with additional information. Defaults to loading
    an spcat catalog.

    Anything in catdict will overwrite what's loaded in from the
    read catalog function, so use cautiously.
    '''

    new_dict = _read_spcat(filein)  # read in the catalog file and produce the
    # dictionary

    cat = Catalog(catdict=new_dict)  # make the catalog and return it

    return cat


def _make_level_dict(qn1low, qn2low, qn3low, qn4low, qn5low,
                     qn6low, qn7low, qn8low, qn1up, qn2up,
                     qn3up, qn4up, qn5up, qn6up, qn7up, qn8up,
                     frequency, elow, gup, qn_list_low, qn_list_up,
                     level_qns, level_dict, qnstrfmt=None):

    # we need to sort out unique levels from our catalog.
    # Those will have unique quantum numbers. When we find
    # a match to a lower level, add the info in.
    for x in range(len(frequency)):
        qnstr_low = qn_list_low[x]
        level_dict[qnstr_low] = {'energy':	elow[x],
                                 'g':	None,
                                 'g_flag':	False,
                                 'qn1':	qn1low[x] if qn1low is not None
                                 else None,
                                 'qn2':	qn2low[x] if qn2low is not None
                                 else None,
                                 'qn3':	qn3low[x] if qn3low is not None
                                 else None,
                                 'qn4':	qn4low[x] if qn4low is not None
                                 else None,
                                 'qn5':	qn5low[x] if qn5low is not None
                                 else None,
                                 'qn6':	qn6low[x] if qn6low is not None
                                 else None,
                                 'qn7':	qn7low[x] if qn7low is not None
                                 else None,
                                 'qn8':	qn8low[x] if qn8low is not None
                                 else None,
                                 'id':	qn_list_low[x],
                                 'qnstrfmt':	qnstrfmt,
                                 }

    # do it again to fill in energy levels that were
    # upper states and didn't get hit
    for x in range(len(frequency)):
        qnstr_up = qn_list_up[x]
        if level_dict[qnstr_up] is None:
            # calculate the energy.  Move the transition from MHz -> cm-1 -> K
            freq_cm = (frequency[x]*1E6/ccm)
            freq_K = freq_cm / kcm
            level_dict[qnstr_up] = {'energy':	elow[x] + freq_K,
                                    'g':	gup[x],
                                    'g_flag':	False,
                                    'qn1':	qn1up[x] if qn1up is not None
                                    else None,
                                    'qn2':	qn2up[x] if qn2up is not None
                                    else None,
                                    'qn3':	qn3up[x] if qn3up is not None
                                    else None,
                                    'qn4':	qn4up[x] if qn4up is not None
                                    else None,
                                    'qn5':	qn5up[x] if qn5up is not None
                                    else None,
                                    'qn6':	qn6up[x] if qn6up is not None
                                    else None,
                                    'qn7':	qn7up[x] if qn7up is not None
                                    else None,
                                    'qn8':	qn8up[x] if qn8up is not None
                                    else None,
                                    'id':	qn_list_up[x],
                                    'qnstrfmt':	qnstrfmt,
                                    }

    # go grab the degeneracies
    for x in range(len(frequency)):
        qnstr_up = qn_list_up[x]
        if level_dict[qnstr_up]['g'] is None:
            level_dict[qnstr_up]['g'] = gup[x]

    # now go through and fill any degeneracies that didn't
    # get hit (probably ground states)
    # assume it's just 2J+1.  Set the flag for a
    # calculated degeneracy to True.
    for x in level_dict:
        if level_dict[x]['g'] is None:
            level_dict[x]['g'] = 2*level_dict[x]['qn1'] + 1
            level_dict[x]['g_flag'] = True

    return level_dict


def parse_cat(filein, qn_label_list, qnstrfmt=None, partition_dict=None,
              qpart_file=None):
    '''
    Loads a molecule in from a catalog file.  Default catalog
    type is molsim.  Override things with catdict.  Generates
    energy level objects, transition objects, a partition
    function object and	a molecule object which it returns.
    '''

    # load the catalog in
    cat = _load_catalog(filein)

    # now we have to make a hash for every entries upper and
    # lower state.  If this already exists in the catalog, great.
    # If not, we have to make it.
    if cat.qnup_str is None:
        qn_list_up = [_make_qnstr(
            qn1, qn2, qn3, qn4, qn5, qn6, qn7, qn8)
            for qn1, qn2, qn3, qn4, qn5, qn6, qn7, qn8 in zip(
            cat.qn1up, cat.qn2up, cat.qn3up, cat.qn4up,
            cat.qn5up, cat.qn6up, cat.qn7up, cat.qn8up)]
    else:
        qn_list_up = cat.qnup_str

    if cat.qnlow_str is None:
        qn_list_low = [_make_qnstr(qn1, qn2, qn3, qn4, qn5, qn6, qn7, qn8)
                       for qn1, qn2, qn3, qn4, qn5, qn6, qn7, qn8 in zip(
            cat.qn1low, cat.qn2low, cat.qn3low, cat.qn4low, cat.qn5low,
            cat.qn6low, cat.qn7low, cat.qn8low)]
    else:
        qn_list_low = cat.qnlow_str

    level_qns = np.concatenate((qn_list_low, qn_list_up), axis=0)
    level_qns = list(set(list(level_qns)))  # get the unique ones
    level_dict = dict.fromkeys(level_qns)

    # now we find unique energy levels.  We just get the dictionary
    # of levels, since that function is computationally intensive
    # and we want to njit it.

    level_dict = _make_level_dict(
        cat.qn1low,
        cat.qn2low,
        cat.qn3low,
        cat.qn4low,
        cat.qn5low,
        cat.qn6low,
        cat.qn7low,
        cat.qn8low,
        cat.qn1up,
        cat.qn2up,
        cat.qn3up,
        cat.qn4up,
        cat.qn5up,
        cat.qn6up,
        cat.qn7up,
        cat.qn8up,
        cat.frequency,
        cat.elow,
        cat.gup,
        qn_list_low,
        qn_list_up,
        level_qns,
        level_dict,
        qnstrfmt,
    )

    # load those levels into level objects and add to a list
    levels = []
    for x in level_dict:
        levels.append(Level(
            energy=level_dict[x]['energy'],
            g=level_dict[x]['g'],
            g_flag=level_dict[x]['g_flag'],
            qn1=level_dict[x]['qn1'],
            qn2=level_dict[x]['qn2'],
            qn3=level_dict[x]['qn3'],
            qn4=level_dict[x]['qn4'],
            qn5=level_dict[x]['qn5'],
            qn6=level_dict[x]['qn6'],
            qn7=level_dict[x]['qn7'],
            qn8=level_dict[x]['qn8'],
            qnstrfmt=level_dict[x]['qnstrfmt'],
            id=level_dict[x]['id']
        ))
    # sort them so the lowest energy is first
    levels.sort(key=lambda x: x.energy)

    # we'll now update the catalog with some of the things we've
    # calculated like eup and glow unless they're already present
    tmp_dict = {}
    for x in levels:
        tmp_dict[x.id] = [x.g, x.energy]
    if cat.glow is None:
        cat.glow = np.empty_like(cat.frequency)
    if cat.eup is None:
        cat.eup = np.empty_like(cat.frequency)
    for x in range(len(cat.frequency)):
        cat.glow[x] = tmp_dict[cat.qnlow_str[x]][0]
        cat.eup[x] = tmp_dict[cat.qnup_str[x]][1]

    # now we have to load the transitions in and make transition objects

    # make a partition function object and assign it to the molecule
    # if there's no other info, assume we're state counting
    if partition_dict is None:
        partition_dict = {}
    # if there's a qpart file specified, read that in
    if qpart_file is not None:
        partition_dict['qpart_file'] = qpart_file
    # make the partition function object and assign it
    qpart = PartitionFunction(
        qpart_file=partition_dict['qpart_file'] if 'qpart_file' in
        partition_dict else None,
        form=partition_dict['form'] if 'form' in partition_dict
        else None,
        params=partition_dict['params'] if 'params' in partition_dict
        else None,
        temps=partition_dict['temps'] if 'temps' in partition_dict
        else None,
        vals=partition_dict['vals'] if 'vals' in partition_dict else None,
        mol=partition_dict['mol'] if 'mol' in partition_dict else None,
        gs=partition_dict['gs'] if 'gs' in partition_dict else None,
        energies=partition_dict['energies'] if 'energies' in partition_dict
        else None,
        sigma=partition_dict['sigma'] if 'sigma' in partition_dict else 1.,
        vib_states=partition_dict['vib_states'] if 'vib_states' in
        partition_dict else None,
        vib_is_K=partition_dict['vib_is_K'] if 'vib_is_K' in partition_dict
        else None,
        notes=partition_dict['notes'] if 'notes' in partition_dict else None,
    )

    # set sijmu and aij
    cat._set_sijmu_aij(qpart)
    if len(qn_label_list) != len(cat.qnlow_str[0])/2 or \
            len(qn_label_list) != len(cat.qnup_str[0])/2:
        raise ValueError(
            'Quantum number labels do not match the number of \
            quantum numbers in .cat file')
    lower_state_qn_dict_list = []
    upper_state_qn_dict_list = []

    # Reformat lower and upper quantum state number strings into dictionaries
    # The dictionaries are appended to lists.
    for i in range(len(cat.frequency)):
        lower_state_qn_dict = {}
        upper_state_qn_dict = {}
        j = 0
        for qn_label in qn_label_list:
            lower_state_qn_dict[qn_label] = int(cat.qnlow_str[i][j:j+2])
            upper_state_qn_dict[qn_label] = int(cat.qnup_str[i][j:j+2])
            j += 2
        lower_state_qn_dict_list.append(lower_state_qn_dict)
        upper_state_qn_dict_list.append(upper_state_qn_dict)

    return cat.frequency, cat.freq_err, cat.logint, cat.sijmu, \
        cat.aij, cat.elow, cat.eup, cat.glow, cat.gup, cat.qnformat, \
        cat.qnlow_str, cat.qnup_str, lower_state_qn_dict_list, \
        upper_state_qn_dict_list
