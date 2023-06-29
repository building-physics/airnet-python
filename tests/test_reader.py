# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import airnet

def test_afdata_short_lines():
    bad_elements = [['element  flow-1.000 cfr', 'element cpt1 plr 0.000399207'],
                    ['node   node-1    c   0.0   20.0', 'element  duct-2x25  dwc  2.0   0.25'],
                    ['element  rest3  plr  1.0e-5  4.092e-06'],
                    ['element  duct-2x25  dwc  2.0   0.25  0.04909'],
                    ['element  duct-2x25  dwc  2.0   0.25  0.04909  0.00015   duct 2m x 25cm'],
                    ['element  duct-3x25  dwc  3.0   0.25  0.04909  0.00015   duct 3m x 25c',
                     '0.0   64.0'],
                    ['element  dor-1.60  dor  0.015575  0.015575   1.76494',
                     '0.0001    2.0       0.8        0.78'],
                    ['element  dor-1.60  dor  0.015575  0.015575   1.76494   .500    doorway - 1.6m^2',
                     '0.0001    2.0       0.8'],
                    ['element  dor-1.60  dor  0.015575  0.015575   1.76494   .500    doorway - 1.6m^2'],
                    ['element  fan-A  cpf  10.0  1.0'],
                    ['element fan-A fan 8.1216e-06 8.1216e-06'],
                    ['element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A',
                     '1.204   764.429   5.46267'],
                    ['element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A',
                     '1.204   764.429   5.46267  0.02   1   10000.'],
                    ['element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A'],
                    ['element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A',
                     '1.204   764.429   5.46267  0.02   1   10000.',
                     '764.429 -18.2932 19.4639'],
                    ['element dqfr-2x25 qfr 0.310242'],
                    ['element  test      ckv  0.0'],
                    ['link   link-1   node-1   0.9   node-2   0.9'],
                    ['link   link-1   node-1   0.0   node-2   0.0   fan-A   notnull']]

    for el in bad_elements:
      reader = airnet.Reader(iter(el))
      results = []

      try:
          for item in reader:
              results.append(item)
      except airnet.BadNetworkInput:
          pass

      assert len(results) == 0

AFDATA_CP1 = '''/*subfile:  afdata.cp1  ******************************************************/
/
title  constant power fan test #1

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0

element  fan-A  cpf  10.0  1.0  2.0
element cpt1 plr 0.000399207 0.000399207 0.756938 0.5 ! 

link   link-1   node-1   0.0   node-2   0.0   fan-A   null
link   link-2   node-2   0.0   node-1   0.0   cpt1    null

*********

'''

def test_afdata_cp1():
    reader = airnet.Reader(iter(AFDATA_CP1.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 7
    assert results[3]['input_type'] == airnet.InputType.ELEMENT
    assert results[6]['input_type'] == airnet.InputType.LINK

AFDATA_CPX = '''/*subfile:  afdata.cp1  ******************************************************/
/
title  constant power fan test #1

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0

element  fan-A  cpx  10.0  1.0  2.0
element cpt1 plr 0.000399207 0.000399207 0.756938 0.5 ! 

link   link-1   node-1   0.0   node-2   0.0   fan-A   null
link   link-2   node-2   0.0   node-1   0.0   cpt1    null

*********

'''

def test_afdata_cpx():
    reader = airnet.Reader(iter(AFDATA_CPX.splitlines()))
    results = []

    try:
        for item in reader:
          results.append(item)
    except airnet.BadNetworkInput:
        pass
    
    assert len(results) == 3
    assert results[2]['input_type'] == airnet.InputType.NODE
    assert results[2]['ht'] == 0.0

AFDATA_CV1 = '''/*subfile:  afdata.cv1  ******************************************************/
/
title  check valve test #1 input file

node   node-1    c   0.0   20.0   16.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0    0.0

element  orf-0.01  plr  7.2e-6  7.2e-6  0.00848528  0.5   orifice - 0.01 m^2
element  test      ckv  0.0     10.0                      check valve

link   link-1   node-1   0.0   node-2   0.0   orf-0.01   null
link   link-2   node-2   0.0   node-3   0.0   test       null

*********

'''

def test_afdata_cv1():
    reader = airnet.Reader(iter(AFDATA_CV1.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 8
    assert results[5]['input_type'] == airnet.InputType.ELEMENT
    assert results[7]['input_type'] == airnet.InputType.LINK

AFDATA_DR1 = '''/*subfile:  afdata.dr1  ******************************************************/
/
title  doorway test #1 input file

node   node-1    c   0.0   18.0    0.0
node   node-2    v   0.0   22.0

element  orf-0.16  plr  0.0015575 0.0015575  0.176494  .500    opening - .16m^2
element  dor-1.60  dor  0.015575  0.015575   1.76494   .500    doorway - 1.6m^2
                        0.0001    2.0       0.8        0.78

link   link-1   node-1   0.9   node-2   0.9   orf-0.16   null
link   link-2   node-1   0.7   node-2   0.7   orf-0.16   null
link   link-3   node-1   0.5   node-2   0.5   orf-0.16   null
link   link-4   node-1   0.3   node-2   0.3   orf-0.16   null
link   link-5   node-1   0.1   node-2   0.1   orf-0.16   null
link   link-6   node-1  -0.1   node-2  -0.1   orf-0.16   null
link   link-7   node-1  -0.3   node-2  -0.3   orf-0.16   null
link   link-8   node-1  -0.5   node-2  -0.5   orf-0.16   null
link   link-9   node-1  -0.7   node-2  -0.7   orf-0.16   null
link   link-10  node-1  -0.9   node-2  -0.9   orf-0.16   null
link   link-11  node-1  -1.0   node-2  -1.0   dor-1.60   null

*********

'''

def test_afdata_dr1():
    reader = airnet.Reader(iter(AFDATA_DR1.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 16
    assert results[2]['input_type'] == airnet.InputType.NODE
    assert results[4]['type'] == airnet.ElementType.DOR

AFDATA_DW1 = '''/*subfile:  afdata.dw1  ******************************************************/
/
title  duct (Darcy-Weisbach-Colebrook model) test #1 input file

node   node-1    c   0.0   20.0    9.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    c   0.0   20.0    0.0

element  duct-2x25  dwc  2.0   0.25  0.04909  0.00015   duct 2m x 25cm
                         0.0   64.0  0.0      128.
element  duct-3x25  dwc  3.0   0.25  0.04909  0.00015   duct 3m x 25cm
                         0.0   64.0  0.0      128.
element  duct-5x25  dwc  5.0   0.25  0.04909  0.00015   duct 5m x 25cm
                         0.0   64.0  0.0      128.
element  duct-10x25 dwc 10.0   0.25  0.04909  0.00015   duct 10m x 25cm
                         0.0   64.0  0.0      128.

link   link-1   node-1   0.0   node-2   0.0   duct-2x25   null
link   link-2   node-2   0.0   node-3   0.0   duct-3x25   null
link   link-3   node-3   0.0   node-4   0.0   duct-5x25   null
link   link-4   node-1   0.0   node-4   0.0   duct-10x25  null

*********

'''

AFDATA_CF1 = '''/*subfile:  afdata.cf1  ******************************************************/
/
title constant flow test #1 input file

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0    0.0

element  orf-0.1000 plr 2.2769e-4 2.2769e-4  0.0848528    0.5  orf - 0.1 m^2
element  flow-1.000 cfr 1.0                       constant flow of 1.000 kg/s

link   link-1   node-1   0.0   node-2   0.0   orf-0.1000   null
link   link-2   node-2   0.0   node-3   0.0   flow-1.000   null

*********

'''

def test_afdata_cf1():
    reader = airnet.Reader(iter(AFDATA_CF1.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 8
    assert results[2]['input_type'] == airnet.InputType.NODE
    assert results[5]['type'] == airnet.ElementType.CFR


AFDATA_FN1 = '''/*subfile:  afdata.fn1  ******************************************************/
/
title  fan test #1 input file:  comparison to Osborne problem #1

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0    0.0

element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A
         1.204   764.429   5.46267  0.02   1   10000.
         764.429 -18.2932 19.4639 -7.63948 -10000.
element  resst  plr  3.0e-6  9.34e-06  0.11199   0.5  ! problem 1 resistance

link   link-1   node-1   0.0   node-2   0.0   fan-A   null
link   link-2   node-2   0.0   node-3   0.0   resst   null

*********


FAN -> FAN:
        Cross sectional area = 0.01 m^2
          Hydraulic diameter = 0.1128 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 0.000303616 Pa
    Laminar flow coefficient = 8.1216e-06 
  Turbulent flow coefficient = 0.00848528 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 0.000134746  0.000162234  0.000162234

         Data points (dP Q F):           765           0           0
         Data points (dP Q F):           755           1       1.204
         Data points (dP Q F):           730           2       2.408
         Data points (dP Q F):           590           3       3.612
         Data points (dP Q F):           275           4       4.816
Polynomial coefficients:  764.429  -18.2932  19.4639  -7.63948
Fit (dP Fdata Ffit diff):           765           0     764.429   -0.571289
Fit (dP Fdata Ffit diff):           755       1.204     757.285      2.2854
Fit (dP Fdata Ffit diff):           730       2.408     726.571    -3.42853
Fit (dP Fdata Ffit diff):           590       3.612     592.286     2.28607
Fit (dP Fdata Ffit diff):           275       4.816     274.428   -0.571594

element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A
               1.204 764.429 5.46267 0.02 1 10000
             764.429 -18.2932 19.4639 -7.63948 -10000

'''

def test_afdata_fn1():
    reader = airnet.Reader(iter(AFDATA_FN1.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 8
    assert results[4]['type'] == airnet.ElementType.FAN
    assert results[5]['expt'] == 0.5

AFDATA_FN2 = '''/*subfile:  afdata.fn2  ******************************************************/
/
title  fan test #2 input file:  comparison to Osborne problem #5

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    v   0.0   20.0

element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A
         1.204   764.429   5.46267  0.02   1   10000.
         764.429 -18.2932 19.4639 -7.63948 -10000.
element  rest1  plr  1.0e-5  1.585e-05  0.19005   0.5   ! resistance 1
element  rest2  plr  1.0e-5  1.110e-05  0.13306   0.5   ! resistance 2
element  rest3  plr  1.0e-5  2.906e-05  0.34843   0.5   ! resistance 3

link   link-1   node-1   0.0   node-2   0.0   fan-A   null
link   link-2   node-1   0.0   node-3   0.0   fan-A   null
link   link-3   node-2   0.0   node-4   0.0   rest1   null
link   link-4   node-3   0.0   node-4   0.0   rest2   null
link   link-5   node-4   0.0   node-1   0.0   rest3   null

*********

'''

AFDATA_FN3 = '''/*subfile:  afdata.fn3  ******************************************************/
/
title  fan test #3 input file:  comparison to Osborne problem #6

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    v   0.0   20.0
node   node-5    c   0.0   20.0    0.0

element fan-A fan 8.1216e-06 8.1216e-06 0.00848528 0.5 ! Osborne fan A
         1.204   764.429   5.46267  0.02   1   10000.
         764.429 -18.2932 19.4639 -7.63948 -10000.
element  rest1  plr  1.0e-5  1.585e-05  0.19005   0.5  ! resistance 1
element  rest2  plr  1.0e-5  1.110e-05  0.13306   0.5  ! resistance 2
element  rest3  plr  1.0e-5  4.092e-06  0.04907   0.5  ! resistance 3

link   link-1   node-1   0.0   node-2   0.0   fan-A   null
link   link-2   node-2   0.0   node-3   0.0   rest1   null
link   link-3   node-3   0.0   node-4   0.0   fan-A   null
link   link-4   node-3   0.0   node-5   0.0   rest3   null
link   link-5   node-4   0.0   node-5   0.0   rest2   null

*********

'''

AFDATA_MR1 = '''/*subfile:  afdata.mr1  ******************************************************/
/
title  multi-range resistance test #1 input file

node   node-1    c   0.0   20.0   500.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    c   0.0   20.0    0.0

element  orf-0.1 plr 0.000256914 0.000256914 0.0848528 0.5 ! orifice - 0.1 m^2
element  mrr1    mrr   5                     -1000         ! test case
                     0.1            0     -1.15356
               0.0931065          0.5 -3.03412e-05
                  16.903            1  3.03412e-05
               0.0931053          0.5       115.36
                       1            0         1000

link   link-1   node-1   0.0   node-2   0.0   orf-0.1    null
link   link-2   node-2   0.0   node-4   0.0   mrr1       null
link   link-3   node-1   0.0   node-3   0.0   orf-0.1    null
link   link-4   node-3   0.0   node-4   0.0   orf-0.1    null

*********


DATA -> MRR:
 #      Q(m^3/s)      F(kg/s)       dP(Pa)        a           b
 0   -0.0830565         -0.1        -1000
 1   -0.0830565         -0.1     -1.15356          0.1            0
 2  -0.00042596 -0.000512856 -3.03412e-05    0.0931065          0.5
 3   0.00042596  0.000512856  3.03412e-05       16.903            1
 4     0.830565            1       115.36    0.0931053     0.499999
 5     0.830565            1         1000            1            0

element test mrr1 5 2 -1000 ! testcase
                     0.1            0     -1.15356
               0.0931065          0.5 -3.03412e-05
                  16.903            1  3.03412e-05
               0.0931053     0.499999       115.36
                       1            0         1000

'''

AFDATA_MR2 = '''/*subfile:  afdata.mr2  ******************************************************/
/
title  multi-range resistance test #2 input file

node   node-1    c   0.0   20.0   10.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0    0.0

element  orf-0.1   plr 2.2769e-4 2.2769e-4  0.0848528    0.5  orf - 0.1 m^2
element  orf-1.0   plr 0.0072    0.0072     0.848528     0.5  orf - 1.0 m^2
element  mrr1      mrr  5          2       -1000              test case
                     0.1           0         -10
               0.0562341        0.25      -0.001
                      10           1       0.001
                0.248322    0.465005          20
                       1           0        1000

link   link-1   node-1   0.0   node-2   0.0   orf-0.1    null
link   link-2   node-2   0.0   node-3   0.0   orf-0.1    null

*********

'''

AFDATA_PL1 = '''/*subfile:  afdata.pl1  ******************************************************/
/
title  powerlaw test #1 input file

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0   -1.0

element  orf-0.01  plr  7.2e-6  7.2e-6  0.00848528  0.5   orifice - 0.01 m^2
element  orf-0.04  plr  5.76e-5 5.76e-5 0.03394113  0.5   orifice - 0.04 m^2

link   link-1   node-1   0.0   node-2   0.0   orf-0.01   null
link   link-2   node-2   0.0   node-3   0.0   orf-0.04   null

*********

'''

def test_afdata_pl1():
    reader = airnet.Reader(iter(AFDATA_PL1.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 8
    assert results[4]['expt'] == 0.5
    assert results[5]['expt'] == 0.5


AFDATA_PL1_TITLE = '''/*subfile:  afdata.pl1  ******************************************************/
/
title  powerlaw test #1 input file

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0   -1.0
title  powerlaw test #1 input file

element  orf-0.01  plr  7.2e-6  7.2e-6  0.00848528  0.5   orifice - 0.01 m^2
element  orf-0.04  plr  5.76e-5 5.76e-5 0.03394113  0.5   orifice - 0.04 m^2

link   link-1   node-1   0.0   node-2   0.0   orf-0.01   null
link   link-2   node-2   0.0   node-3   0.0   orf-0.04   null

*********

'''

def test_afdata_pl1_title():
    reader = airnet.Reader(iter(AFDATA_PL1_TITLE.splitlines()))
    results = []

    try:
      for item in reader:
          results.append(item)
    except airnet.BadNetworkInput:
        pass

    assert len(results) == 4
    assert results[0]['input_type'] == airnet.InputType.TITLE
    assert results[1]['input_type'] == airnet.InputType.NODE

AFDATA_PL1_BAD_NODE_1 = '''/*subfile:  afdata.pl1  ******************************************************/
/
title  powerlaw test #1 input file

node   node-1    c   0.0   20.0    0.0
node   node-2    x   0.0   20.0
node   node-3    c   0.0   20.0   -1.0
title  powerlaw test #1 input file

element  orf-0.01  plr  7.2e-6  7.2e-6  0.00848528  0.5   orifice - 0.01 m^2
element  orf-0.04  plr  5.76e-5 5.76e-5 0.03394113  0.5   orifice - 0.04 m^2

link   link-1   node-1   0.0   node-2   0.0   orf-0.01   null
link   link-2   node-2   0.0   node-3   0.0   orf-0.04   null

*********

'''

def test_afdata_pl1_bad_node_1():
    reader = airnet.Reader(iter(AFDATA_PL1_BAD_NODE_1.splitlines()))
    results = []

    try:
      for item in reader:
          results.append(item)
    except airnet.BadNetworkInput:
        pass

    assert len(results) == 2
    assert results[0]['input_type'] == airnet.InputType.TITLE
    assert results[1]['input_type'] == airnet.InputType.NODE

AFDATA_PL1_BAD_NODE_2 = '''/*subfile:  afdata.pl1  ******************************************************/
/
title  powerlaw test #1 input file

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0
node   node-3    c   0.0   20.0   -1.0
title  powerlaw test #1 input file

element  orf-0.01  plr  7.2e-6  7.2e-6  0.00848528  0.5   orifice - 0.01 m^2
element  orf-0.04  plr  5.76e-5 5.76e-5 0.03394113  0.5   orifice - 0.04 m^2

link   link-1   node-1   0.0   node-2   0.0   orf-0.01   null
link   link-2   node-2   0.0   node-3   0.0   orf-0.04   null

*********

'''

def test_afdata_pl1_bad_node_2():
    reader = airnet.Reader(iter(AFDATA_PL1_BAD_NODE_2.splitlines()))
    results = []

    try:
      for item in reader:
          results.append(item)
    except airnet.BadNetworkInput:
        pass

    assert len(results) == 2
    assert results[0]['input_type'] == airnet.InputType.TITLE
    assert results[1]['input_type'] == airnet.InputType.NODE

AFDATA_PL2 = '''/*subfile:  afdata.pl2  ******************************************************/
/
title  powerlaw test #2 input file

node   node-1    c   0.0   20.0   0.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    c   0.0   20.0  -100.0

element  orf-0.0001 plr 8.124e-09 8.124e-09 8.48528e-05 0.5 ! orf - 0.0001 m^2
element  orf-0.001  plr 2.569e-07 2.569e-07 0.000848528 0.5 ! orf - 0.001 m^2
element  orf-0.01   plr 8.124e-06 8.124e-06 0.00848528  0.5 ! orf - 0.01 m^2
element  orf-0.1    plr 0.0002569 0.0002569 0.0848528   0.5 ! orf - 0.1 m^2
element  orf-1.0    plr 0.008124  0.008124  0.848528    0.5 ! orf - 1.0 m^2
element  orf-10.0   plr 0.2569    0.2569    8.48528     0.5 ! orf - 10.0 m^2
element  orf-100.0  plr 8.124     8.124     84.8528     0.5 ! orf - 100.0 m^2

link   link-1   node-1   0.0   node-2   0.0   orf-0.0001   null
link   link-2   node-2   0.0   node-3   0.0   orf-0.0001   null
link   link-3   node-3   0.0   node-4   0.0   orf-0.0001   null

*********

ORIFICE -> PLR:
        Cross sectional area = 0.0001 m^2
          Hydraulic diameter = 0.0112838 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 0.0303412 Pa
    Laminar flow coefficient = 8.12434e-09 
  Turbulent flow coefficient = 8.48528e-05 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 1.347e-05  1.62179e-05  1.62179e-05

element orf-0.0001 plr 8.12434e-09 8.12434e-09 8.48528e-05 0.5 ! orifice - 0.0001 m^2


ORIFICE -> PLR:
        Cross sectional area = 0.001 m^2
          Hydraulic diameter = 0.0356825 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 0.00303412 Pa
    Laminar flow coefficient = 2.56914e-07 
  Turbulent flow coefficient = 0.000848528 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 4.2596e-05  5.12856e-05  5.12856e-05

element orf-0.001 plr 2.56914e-07 2.56914e-07 0.000848528 0.5 ! orifice - 0.001 m^2


ORIFICE -> PLR:
        Cross sectional area = 0.01 m^2
          Hydraulic diameter = 0.112838 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 0.000303412 Pa
    Laminar flow coefficient = 8.12434e-06 
  Turbulent flow coefficient = 0.00848528 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 0.000134701  0.000162179  0.000162179

element orf-0.01 plr 8.12434e-06 8.12434e-06 0.00848528 0.5 ! orifice - 0.01 m^2


ORIFICE -> PLR:
        Cross sectional area = 0.1 m^2
          Hydraulic diameter = 0.356825 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 3.03412e-05 Pa
    Laminar flow coefficient = 0.000256914 
  Turbulent flow coefficient = 0.0848528 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 0.00042596  0.000512856  0.000512856

element orf-0.1 plr 0.000256914 0.000256914 0.0848528 0.5 ! orifice - 0.1 m^2


ORIFICE -> PLR:
        Cross sectional area = 1 m^2
          Hydraulic diameter = 1.12838 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 3.03412e-06 Pa
    Laminar flow coefficient = 0.00812434 
  Turbulent flow coefficient = 0.848528 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 0.00134701  0.00162179  0.00162179

element orf-1.0 plr 0.00812434 0.00812434 0.848528 0.5 ! orifice - 1.0 m^2


ORIFICE -> PLR:
        Cross sectional area = 10 m^2
          Hydraulic diameter = 3.56825 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 3.03412e-07 Pa
    Laminar flow coefficient = 0.256914 
  Turbulent flow coefficient = 8.48528 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 0.0042596  0.00512856  0.00512856

element orf-10.0 plr 0.256914 0.256914 8.48528 0.5 ! orifice - 10.0 m^2


ORIFICE -> PLR:
        Cross sectional area = 100 m^2
          Hydraulic diameter = 11.2838 m
       Discharge coefficient = 0.6
  Transition Reynolds number = 100
Laminar/turbulent transition = 3.03412e-08 Pa
    Laminar flow coefficient = 8.12434 
  Turbulent flow coefficient = 84.8528 
Pressure difference exponent = 0.5 
Transition flows (Qre, Fl, Ft): 0.01347  0.0162179  0.0162179

element orf-100.0 plr 8.12434 8.12434 84.8528 0.5 ! orifice - 100.0 m^2


'''

def test_afdata_pl2():
    reader = airnet.Reader(iter(AFDATA_PL2.splitlines()), floats=False)
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 15
    assert results[5]['expt'] == '0.5'
    assert results[6]['expt'] == '0.5'

AFDATA_PL3 = '''/*subfile:  afdata.pl3  ******************************************************/
/
title  powerlaw test #3 input file

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    v   0.0   20.0
node   node-5    v   0.0   20.0
node   node-6    v   0.0   20.0
node   node-7    v   0.0   20.0
node   node-8    v   0.0   20.0
node   node-9    v   0.0   20.0
node   node-10   v   0.0   20.0
node   node-11   v   0.0   20.0
node   node-12   c   0.0   20.0  -100.0

element  orf-0.0001 plr 7.2e-9    7.2e-9     8.48528e-5   0.5  orf - 0.0001 m^2
element  orf-0.0004 plr 5.76e-8   5.76e-8    3.39411e-4   0.5  orf - 0.0004 m^2
element  orf-0.0009 plr 1.944e-7  1.944e-7   7.63675e-4   0.5  orf - 0.0009 m^2
element  orf-0.0025 plr 9.0e-7    9.0e-7     0.00212132   0.5  orf - 0.0025 m^2  
element  orf-0.0036 plr 1.555e-6  1.555e-6   0.00305470   0.5  orf - 0.0036 m^2  
element  orf-0.0100 plr 7.2e-6    7.2e-6     0.00848528   0.5  orf - 0.01 m^2
element  orf-1.0000 plr 0.0072    0.0072     0.848528     0.5  orf - 1.0 m^2
element  orf-4.0000 plr 0.0576    0.0576     3.394113     0.5  orf - 4.0 m^2

link  link-1  node-2  0.0  node-3  0.0  orf-0.0100  null
link  link-2  node-2  0.0  node-3  0.0  orf-1.0000  null
link  link-3  node-2  0.0  node-3  0.0  orf-4.0000  null  
link  link-4  node-3  0.0  node-4  0.0  orf-0.0025  null
link  link-5  node-3  0.0  node-4  0.0  orf-0.0036  null  
link  link-6  node-4  0.0  node-11 0.0  orf-1.0000  null
link  link-7  node-4  0.0  node-11 0.0  orf-1.0000  null  
link  link-8  node-4  0.0  node-11 0.0  orf-1.0000  null
link  link-9  node-1  0.0  node-2  0.0  orf-1.0000  null  
link  link-10 node-2  0.0  node-5  0.0  orf-1.0000  null  
link  link-11 node-5  0.0  node-6  0.0  orf-0.0001  null
link  link-12 node-6  0.0  node-7  0.0  orf-1.0000  null  
link  link-13 node-7  0.0  node-8  0.0  orf-0.0004  null
link  link-14 node-8  0.0  node-11 0.0  orf-4.0000  null  
link  link-15 node-2  0.0  node-9  0.0  orf-0.0004  null
link  link-16 node-9  0.0  node-10 0.0  orf-4.0000  null  
link  link-17 node-9  0.0  node-10 0.0  orf-0.0001  null  
link  link-18 node-9  0.0  node-10 0.0  orf-1.0000  null
link  link-19 node-10 0.0  node-11 0.0  orf-0.0009  null  
link  link-20 node-11 0.0  node-12 0.0  orf-1.0000  null  

*********

'''

AFDATA_PR1 = '''/*subfile:  afdata.pr1  ******************************************************/
/
title  pressure relief valve test #1 input file

node   node-1    c   0.0   20.0   10.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0    0.0

element  orf-0.1   plr 2.2769e-4 2.2769e-4  0.0848528    0.5  orf - 0.1 m^2
element  orf-1.0   plr 0.0072    0.0072     0.848528     0.5  orf - 1.0 m^2
element  test      prv 1.0       10.        0.1               relief valve

link   link-1   node-1   0.0   node-2   0.0   orf-1.0    null
link   link-2   node-2   0.0   node-3   0.0   test       null

*********

'''

AFDATA_PR2 = '''/*subfile:  afdata.pr2  ******************************************************/
/
title  pressure relief valve test #2 input file

node   node-1    c   0.0   20.0  100.0
node   node-2    v   0.0   20.0
node   node-3    c   0.0   20.0    0.0

element  orf-0.1   plr 2.2769e-4 2.2769e-4  0.0848528    0.5  orf - 0.1 m^2
element  orf-1.0   plr 0.0072    0.0072     0.848528     0.5  orf - 1.0 m^2
element  test      prv 1.0       1.0        0.1               relief valve

link   link-1   node-1   0.0   node-2   0.0   orf-0.1    null
link   link-2   node-2   0.0   node-3   0.0   test       null

*********

'''

AFDATA_QF1 = '''/*subfile:  afdata.qf1  ******************************************************/
/
title  quadratic flow test #1 input file - QFR, PLR, & DWC duct models

node   source    c   0.0   20.0    0.01
node   duct-1    v   0.0   20.0
node   duct-2    v   0.0   20.0
node   dplr-1    v   0.0   20.0
node   dplr-2    v   0.0   20.0
node   dqfr-1    v   0.0   20.0
node   dqfr-2    v   0.0   20.0
node   sink      c   0.0   20.0    0.0

element duct-2x25 dwc 2 0.25 0.0490874 0.0015 duct - 2m long by 25cm dia
                      0 64   0 64
element duct-3x25 dwc 3 0.25 0.0490874 0.0015 duct - 3m long by 25cd dia
                      0 64   0 64
element duct-5x25 dwc 5 0.25 0.0490874 0.0015 duct - 5m long by 25cm dia
                      0 64   0 64
element duct-10x25 dwc 10 0.25 0.0490874 0.0015 duct - 10m long by 25cm dia
                      0 64   0 64
element dplr-2x25 plr 3.40642e-05 3.40642e-05 0.133079 0.524431 plr-duct - 2m long by 25cm dia
element dplr-3x25 plr 2.27095e-05 2.27095e-05 0.107587 0.524431 plr-duct - 3m long by 25cm dia
element dplr-5x25 plr 1.36257e-05 1.36257e-05 0.0823032 0.524431 plr-duct - 5m long by 25cm dia
element dqfr-2x25 qfr 0.310242 44.8089 qfr-duct - 2m long by 25cm dia
element dqfr-3x25 qfr 0.465364 67.2134 qfr-duct - 3m long by 25cm dia
element dqfr-5x25 qfr 0.775606 112.022 qfr-duct - 5m long by 25cm dia

link   duct-0   source   0.0   sink     0.0   duct-10x25  null
link   duct-1   source   0.0   duct-1   0.0   duct-2x25   null
link   duct-2   duct-1   0.0   duct-2   0.0   duct-3x25   null
link   duct-3   duct-2   0.0   sink     0.0   duct-5x25   null
link   dplr-1   source   0.0   dplr-1   0.0   dplr-2x25   null
link   dplr-2   dplr-1   0.0   dplr-2   0.0   dplr-3x25   null
link   dplr-3   dplr-2   0.0   sink     0.0   dplr-5x25   null
link   dqfr-1   source   0.0   dqfr-1   0.0   dqfr-2x25   null
link   dqfr-2   dqfr-1   0.0   dqfr-2   0.0   dqfr-3x25   null
link   dqfr-3   dqfr-2   0.0   sink     0.0   dqfr-5x25   null

*********

'''

def test_afdata_qf1():
    reader = airnet.Reader(iter(AFDATA_QF1.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 29

AFDATA_TST = '''/*subfile:  afdata.tst  ******************************************************/
/
all elements input/output test file

node   node-1    c   0.0   20.0    0.0   fan inlet
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.5    warm room
node   node-4    v   0.0   19.5    cool room
node   node-5    v   0.0   20.0
node   node-6    v   0.0   20.0
node   node-7    c   0.0   20.0    0.0   constant flow outlet
node   node-8    v   0.0   20.0
node   node-9    c   0.0   20.0    0.0

element  orf-0.0001 plr 7.2e-9    7.2e-9     8.48528e-5   0.5  orf - 0.0001 m^2
element  orf-0.0004 plr 5.76e-8   5.76e-8    3.39411e-4   0.5  orf - 0.0004 m^2
element  orf-4.0000 plr 0.0576    0.0576     3.394113     0.5  orf - 4.0 m^2
element  fan-A  fan  3.0e-5    7.2e-6  0.084853  .500    Osborne fan A
         1.204  750.0  5.00  0.10   3  -100.0
         764.429  -18.2922  19.4633  -7.63940  1.0       no contraflexure
         764.429  -18.2922  19.4633  -7.63940  5.0
         764.429  -18.2922  19.4633  -7.63940  100.
element  fan-B  cpf  100.                                constant power fan
element  dor-1.60  dor  0.015575  0.015575   1.76494   .500    doorway - 1.6m^2
                        0.0001    2.0       0.8        0.78
element duct-5x25 dwc 5 0.25 0.0490874 0.0015 duct - 5m long by 25cm dia
                      0 64   0 64
element  flow-1.000 cfr  1.0                      constant flow of 1.000 kg/s
element  flow-10.00 cfr 10.0                      constant flow of 10.00 kg/s
element  dqfr-5x25 qfr 0.775606 112.022 qfr-duct - 5m long by 25cm dia
element  check     ckv  0.0     10.0                      check valve
element  mrr1      mrr  5          2    -1000.                test case
                     0.1           0.0     -1.15356
                     0.0931064     0.5     -3.86294e-5
                    14.9803        1.0     +3.86294e-5
                     0.0931064     0.5    115.356
                     1.0           0.0   1000.

link  link-1  node-1  0.0  node-2  0.0  fan-A       null
link  link-2  node-2  0.0  node-3  0.0  orf-0.0001  null
link  link-3  node-3  0.0  node-4  0.0  dor-1.60    null  
link  link-4  node-4  0.0  node-5  0.0  orf-0.0004  null
link  link-5  node-2  0.0  node-5  0.0  duct-5x25   null  
link  link-6  node-6  0.0  node-2  0.0  check       null
link  link-7  node-5  0.0  node-6  0.0  check       null  
link  link-8  node-5  0.0  node-6  0.0  orf-0.0004  null  
link  link-9  node-6  0.0  node-7  0.0  flow-1.000  null
link  link-10 node-5  0.0  node-8  0.0  dqfr-5x25   null  
link  link-11 node-8  0.0  node-9  0.0  mrr1        null  
link  link-12 node-8  0.0  node-9  0.0  orf-0.0004  null

*********

all elements input/output test file

node   node-1    c   0.0   20.0    0.0   fan inlet
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.5    warm room
node   node-4    v   0.0   19.5    cool room
node   node-5    v   0.0   20.0
node   node-6    v   0.0   20.0
node   node-7    c   0.0   20.0    0.0   constant flow outlet
node   node-8    v   0.0   20.0
node   node-9    c   0.0   20.0    0.0

element  orf-0.0001 plr 7.2e-9    7.2e-9     8.48528e-5   0.5  orf - 0.0001 m^2
element  orf-0.0004 plr 5.76e-8   5.76e-8    3.39411e-4   0.5  orf - 0.0004 m^2
element  orf-4.0000 plr 0.0576    0.0576     3.394113     0.5  orf - 4.0 m^2
element  fan-A  fan  3.0e-5    7.2e-6  0.084853  .500    Osborne fan A
         1.204  750.0  5.00  0.10   3  -100.0
         764.429  -18.2922  19.4633  -7.63940  1.0       no contraflexure
         764.429  -18.2922  19.4633  -7.63940  5.0
         764.429  -18.2922  19.4633  -7.63940  100.
element  fan-B  cpf  100.                                constant power fan
element  dor-1.60  dor  0.015575  0.015575   1.76494   .500    doorway - 1.6m^2
                        0.0001    2.0       0.8        0.78
element duct-5x25 dwc 5 0.25 0.0490874 0.0015 duct - 5m long by 25cm dia
                      0 64   0 64
element  flow-1.000 cfr 1.0                       constant flow of 1.000 kg/s
element  flow-10.00 cfr 10.0                      constant flow of 10.00 kg/s
element  dqfr-5x25 qfr 0.775606 112.022 qfr-duct - 5m long by 25cm dia
element  check     ckv  0.0     10.0                      check valve
element  mrr1      mrr  5          2    -1000.                test case
                     0.1           0.0     -1.15356
                     0.0931064     0.5     -3.86294e-5
                    14.9803        1.0     +3.86294e-5
                     0.0931064     0.5    115.356
                     1.0           0.0   1000.

link  link-1  node-1  0.0  node-2  0.0  flow-10.00  null
link  link-2  node-2  0.0  node-3  0.0  orf-0.0001  null
link  link-3  node-3  0.0  node-4  0.0  dor-1.60    null  
link  link-4  node-4  0.0  node-5  0.0  orf-0.0004  null
link  link-5  node-2  0.0  node-5  0.0  duct-5x25   null  
link  link-6  node-6  0.0  node-2  0.0  check       null
link  link-7  node-5  0.0  node-6  0.0  check       null  
link  link-8  node-5  0.0  node-6  0.0  orf-0.0004  null  
link  link-9  node-6  0.0  node-7  0.0  flow-1.000  null
link  link-10 node-5  0.0  node-8  0.0  dqfr-5x25   null  
link  link-11 node-8  0.0  node-9  0.0  mrr1        null  
link  link-12 node-8  0.0  node-9  0.0  orf-0.0004  null

*********

'''

AFDATA_WP2 = '''/*subfile:  afdata.wp2  ******************************************************/
/
title  wind pressure test #2 input file

node   node-1    c   0.0   20.0    0.0
node   node-2    v   0.0   20.0
node   node-3    v   0.0   20.0
node   node-4    c   0.0   20.0    0.0

element  orf-0.001  plr 2.2769e-7 2.2769e-7  0.000848528  0.5  orf - 0.001 m^2
element  orf-0.010  plr 7.2e-6    7.2e-6     0.00848528   0.5  orf - 0.01 m^2

link   link-1   node-1   0.0   node-2   0.0   orf-0.001   north   2.0
link   link-2   node-2   0.0   node-3   0.0   orf-0.010   null
link   link-3   node-4   0.0   node-3   0.0   orf-0.001   south   2.0

*********

'''
def test_afdata_wp2():
    reader = airnet.Reader(iter(AFDATA_WP2.splitlines()))
    results = []

    for item in reader:
        results.append(item)

    assert len(results) == 10
    assert results[5]['expt'] == 0.5
    assert results[6]['expt'] == 0.5
