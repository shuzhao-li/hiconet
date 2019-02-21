Directories
===========
/datasrc - unmodified data retrieved from ImmPort or elsewhere
/archived - files not to be used, intermediate or older versions
/unstructured - files not used in the current automated mining, but possible useful for other purposes

This example contains transcriptomics and antibody data, as listed in project.yaml.


Data format/files for each data type:

  file_data_matrix: 'SDY80_pH1N1_2009_Transcriptomics.txt'
  file_feature_annotation: '' [Optional]
  file_observation_annotation: 'biosample.txt'

All above files use the 1st column as primary ID. 
The 1st row in the data matrix for observation ID.
The 1st column in the data matrix for feature ID; data start 2nd row, 2nd column.
All extra annotation and meta data go into the two annotation files.
Annotation files can be shared/reused between different data types.

Transcriptomics data is good for data matrix
============================================
Watch out for length of header (R convention on one less cell for header).

Transcriptomics data - retrieved via ImmuneSpaceR
SDY80_pH1N1_2009_Transcriptomics.txt

# ImmPort sample to subject mapping:
'SDY80-DR25_Tab/Tab/biosample.txt'

Note:
the goal is to match subjects and time points. It’s often that not the same samples were used for different assays.


What’s in biosample.txt -

BIOSAMPLE_ACCESSION	DESCRIPTION	NAME	PLANNED_VISIT_ACCESSION	STUDY_ACCESSION	STUDY_TIME_COLLECTED	STUDY_TIME_COLLECTED_UNIT	STUDY_TIME_T0_EVENT	STUDY_TIME_T0_EVENT_SPECIFY	SUBJECT_ACCESSION	SUBTYPE	TYPE	WORKSPACE_IDBS630798	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114454		RNA	2771BS630803	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114459		RNA	2771BS630806	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114462		RNA	2771BS630813	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114469		RNA	2771BS630815	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114471		RNA	2771BS630816	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114472		RNA	2771BS630822	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114478		RNA	2771BS630833	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114489		RNA	2771BS630836	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114492		RNA	2771BS630841	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114497		RNA	2771BS630847	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114503		RNA	2771BS630849	Before vaccine		PV1834	SDY80	0	Days	Time of initial vaccine administration		SUB114505		RNA	2771BS630858	After vaccine		PV1835	SDY80	1	Days	Time of initial vaccine administration		SUB114451		RNA	2771

In [76]: # create dict biosample - (subj, day)

In [77]: BS2subj = {}

In [78]: for line in open(sample_sub).readlines()[1:]:
    ...:     a = line.split('\t')
    ...:     BS2subj[a[0]] = (a[9], a[5])
    ...:     

In [79]: BS2subj.items()[:5]
Out[79]: 
[('BS808328', ('SUB114456', '70')),
 ('BS808329', ('SUB114457', '70')),
 ('BS630828', ('SUB114484', '0')),
 ('BS630829', ('SUB114485', '0')),
 ('BS630826', ('SUB114482', '0'))] 

Note:
`biosample.txt` can be used as file_observation_annotation.




Processing antibody data
========================

➜  example_SDY80 ls
README.txt    archived      biosample.txt datasrc       project.yaml  unstructured
➜  example_SDY80 ipython
Python 2.7.13 |Anaconda custom (x86_64)| (default, Dec 20 2016, 23:05:08) 
Type "copyright", "credits" or "license" for more information.

IPython 5.1.0 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.

In [1]: # cleaning and reformatting antibody titer data

In [2]: abfile = 'datasrc/neut_ab_titer_2018-04-02_07-00-43.tsv'

In [3]: w = open(abfile).readlines()

In [4]: virus = [x.split('\t')[7] for x in w[1:]]

In [5]: set(virus)
Out[5]: 
{'A/Brisbane/59/2007',
 'A/California/7/2009',
 'A/Uruguay/716/2007',
 'A_Ca_07_swine',
 'A_Uruguay_716_2007',
 'B/Brisbane/60/2008',
 'B_Brisbane_60_2001'}

In [6]: # need standardize viral strains

In [7]: %paste
strain_dict = dict( zip(['A/Brisbane/59/2007',
    ...:  'A/California/7/2009',
    ...:  'A/Uruguay/716/2007',
    ...:  'A_Ca_07_swine',
    ...:  'A_Uruguay_716_2007',
    ...:  'B/Brisbane/60/2008',
    ...:  'B_Brisbane_60_2001'], ['Brisbane_A', 'California', 'Uruguay', 'California', 'Uruguay', 'Brisbane_B', 'Brisbane_B'] ))

## -- End pasted text --

In [8]: w[2]
Out[8]: 'SUB114467.80\t33.0\tMale\tWhite\tCohort2\t0.0\tDays\tA_Uruguay_716_2007\t320.0\t\n'

In [9]: w[2].split('\t')
Out[9]: 
['SUB114467.80',
 '33.0',
 'Male',
 'White',
 'Cohort2',
 '0.0',
 'Days',
 'A_Uruguay_716_2007',
 '320.0',
 '\n']

In [10]: # there's no real sample ID, only subject and time point

In [11]: # need feature x observation as data matrix; and annotation table for observations,

In [12]: # use subject_days as observation ID

In [13]: # do annotation table

In [14]: observation_ann = []

In [15]: for line in w[1:]:
    ...:     a = line.split('\t')
    ...:     # a[0] is subj; a[5] timepoint
    ...:     observation_ann.append((a[0].replace('.80', '')+'_'+a[5], line))
    ...:     

In [16]: w[0]
Out[16]: 'Participant ID\tAge Reported\tGender\tRace\tCohort\tStudy Time Collected\tStudy Time Collected Unit\tVirus\tValue Reported\tUnit Reported\n'

In [17]: s = 'observation_ID\t' + '\t'.join(w[0].split('\t')[:6])

In [18]: s
Out[18]: 'observation_ID\tParticipant ID\tAge Reported\tGender\tRace\tCohort\tStudy Time Collected'

In [19]: s = 'observation_ID\t' + '\t'.join(w[0].split('\t')[:7])

In [20]: s
Out[20]: 'observation_ID\tParticipant ID\tAge Reported\tGender\tRace\tCohort\tStudy Time Collected\tStudy Time Collected Unit'

In [21]: observation_ann2 = []

In [22]: for L in observation_ann:
    ...:     observation_ann2.append( '\t'.join([L[0]] + L[1].split('\t')[:7]) )
    ...:     

In [23]: len(observation_ann2), len(set(observation_ann2))
Out[23]: (960, 240)

In [24]: # so there's 240 observations (4 features each, i.e. 4 viral strains)

In [25]: s += '\n'

In [26]: for L in set(observation_ann2): s += L + '\n'

In [27]: with open('neut_ab_titer__observation_annotation.txt', 'w') as O:
    ...:     O.write(s)
    ...:     

In [28]: 

In [28]: # construct data matrix from observation_ann

In [29]: wanted_strains = ['Brisbane_A', 'Brisbane_B', 'California', 'Uruguay']

In [30]: my_dict = {}

In [31]: for L in observation_ann:
    ...:     a = L[1].split('\t')
    ...:     my_dict[(L[0], strain_dict[a[7]])] = a[8]
    ...:     
    ...:     

In [32]: len(my_dict)
Out[32]: 960

In [33]: # write out data matrix as features x observations

In [34]: wanted_obs = list(set([L[0] for L in observation_ann]))

In [35]: len(wanted_obs)
Out[35]: 240

In [36]: wanted_obs[:3]
Out[36]: ['SUB114500_7.0', 'SUB114499_-7.0', 'SUB114486_0.0']

In [37]: s = 'neut_ab_titer\t' + '\t'.join(wanted_obs) + '\n'

In [38]: for i in wanted_strains:
    ...:     tmp = [my_dict.get((j,i), '') for j in wanted_obs]
    ...:     s += '\t'.join( [i] + tmp ) + '\n'
    ...:     

In [39]: with open('neut_ab_titer__data_matrix.txt', 'w') as O:
    ...:     O.write(s)
    ...:     

In [40]: # inspect written files in Spreadsheets


