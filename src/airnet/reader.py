# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import enum

class BadNetworkInput(Exception):
    pass

class InputType(enum.Enum):
    TITLE = 0
    NODE = 1
    ELEMENT = 2
    LINK = 3

class ElementType(enum.Enum):
    PLR = 0
    DWC = 1
    DOR = 2
    CFR = 3
    FAN = 4
    CPF = 5
    QFR = 6
    CKV = 7

def noop(input):
    return input

class Reader:
    def __init__(self, fp, floats=True):
        self.fp = fp
        self.title = None
        self.line_number = 0
        self.handle_float = noop
        if floats:
            self.handle_float = float

    def get_next(self):
        next_line = next(self.fp)
        self.line_number += 1
        while(next_line.strip() == ''):
            next_line = next(self.fp)
            self.line_number += 1
        return next_line
    
    def __iter__(self):
        return self

    def __next__(self):
        next_line = self.get_next()
        while not next_line.startswith('*'):
            if next_line.lstrip().startswith('title'):
                if self.title:
                    raise BadNetworkInput('Found additional title at line %d' % self.line_number)
                else:
                    self.title = next_line.lstrip().replace('title', '').lstrip()
                    return {'input_type': InputType.TITLE, 'title': self.title}
            data = next_line.split()
            if data[0] == 'node':
                # node name type ht temp pres
                if not data[2] in ['v', 'c', 'a']:
                    raise BadNetworkInput('Node type "%s" at line %d is unrecognized, must be "v", "c", or "a"' % (data[2], self.line_number))
                if data[2] == 'v':
                    if len(data) < 5:
                        raise BadNetworkInput('Node at line %d has fewer than 5 fields' % self.line_number)
                    return {'input_type': InputType.NODE, 'name': data[1], 'type': data[2], 'ht': self.handle_float(data[3]), 
                            'temp': self.handle_float(data[4])}
                else:
                    if len(data) < 6:
                        raise BadNetworkInput('Node at line %d has fewer than 6 fields' % self.line_number)
                    return {'input_type': InputType.NODE, 'name': data[1], 'type': data[2], 'ht': self.handle_float(data[3]), 
                            'temp': self.handle_float(data[4]), 'pres': self.handle_float(data[5])}
            elif data[0] == 'element':
                # element name plr init lam turb expt
                if data[2] == 'plr':
                    if len(data) < 7:
                        raise BadNetworkInput('Element type "plr" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    return {'input_type': InputType.ELEMENT, 'type': ElementType.PLR, 'name': data[1], 'init': self.handle_float(data[3]), 
                            'lam': self.handle_float(data[4]), 'turb': self.handle_float(data[5]), 'expt': self.handle_float(data[6])}
                # element name dwc len dh area rgh
                #         tdlc lflc ldlc init
                elif data[2] == 'dwc':
                    if len(data) < 7:
                        raise BadNetworkInput('Element type "dwc" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    obj = {'input_type': InputType.ELEMENT, 'type': ElementType.DWC, 'name': data[1], 'len': self.handle_float(data[3]), 
                        'dh': self.handle_float(data[4]), 'area': self.handle_float(data[5]), 'rgh': self.handle_float(data[6])}
                    try:
                        next_line = self.get_next()
                    except StopIteration:
                        raise BadNetworkInput('Element type "dwc" at line %d has only one line and cannot be a legal element' % self.line_number)
                    data = next_line.split()
                    if len(data) < 4:
                        raise BadNetworkInput('Element type "dwc" with second line at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    obj['tdlc'] = self.handle_float(data[0]) 
                    obj['lflc'] = self.handle_float(data[1])
                    obj['ldlc'] = self.handle_float(data[2]) 
                    obj['init'] = self.handle_float(data[3])
                    return obj
                # element name dor init lam turb expt
                #         dtmin ht wd cd
                elif data[2] == 'dor':
                    if len(data) < 7:
                        raise BadNetworkInput('Element type "dor" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    obj = {'input_type': InputType.ELEMENT, 'type': ElementType.DOR, 'name': data[1], 'init': self.handle_float(data[3]), 
                        'lam': self.handle_float(data[4]), 'turb': self.handle_float(data[5]), 'expt': self.handle_float(data[6])}
                    try:
                        next_line = self.get_next()
                    except StopIteration:
                        raise BadNetworkInput('Element type "dor" at line %d has only one line and cannot be a legal element' % self.line_number)
                    data = next_line.split()
                    if len(data) < 4:
                        raise BadNetworkInput('Element type "dor" with second line at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    obj['dtmin'] = self.handle_float(data[0]) 
                    obj['ht'] = self.handle_float(data[1])
                    obj['wd'] = self.handle_float(data[2]) 
                    obj['cd'] = self.handle_float(data[3])
                    return obj
                # element name cfr flow
                elif data[2] == 'cfr':
                    if len(data) < 4:
                        raise BadNetworkInput('Element type "cfr" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    return {'input_type': InputType.ELEMENT, 'type': ElementType.CFR, 'name': data[1], 'flow': self.handle_float(data[3])}
                # element name fan init lam turb expt
                #         rdens fdf sop ltt nr mfl
                #         all al2 al3 a14 mf2
                #         a2l a22 a23 a24 mf3
                #         ... ... ... ... mfn
                elif data[2] == 'fan':
                    if len(data) < 7:
                        raise BadNetworkInput('Element type "fan" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    obj = {'input_type': InputType.ELEMENT, 'type': ElementType.FAN, 'name': data[1], 'init': self.handle_float(data[3]), 
                        'lam': self.handle_float(data[4]), 'turb': self.handle_float(data[5]), 'expt': self.handle_float(data[6])}
                    try:
                        next_line = self.get_next()
                    except StopIteration:
                        raise BadNetworkInput('Element type "fan" at line %d has only one line and cannot be a legal element' % self.line_number)
                    data = next_line.split()
                    if len(data) < 6:
                        raise BadNetworkInput('Element type "dwc" with second line at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    obj['rdens'] = self.handle_float(data[0]) 
                    obj['fdf'] = self.handle_float(data[1])
                    obj['sop'] = self.handle_float(data[2]) 
                    obj['ltt'] = self.handle_float(data[3])
                    nr = int(data[4]) 
                    obj['mfl'] = self.handle_float(data[5])
                    pts = []
                    for i in range(nr):
                        try:
                            next_line = self.get_next()
                        except StopIteration:
                            raise BadNetworkInput('Element type "fan" at line %d has insufficient lines of data and cannot be a legal element' % self.line_number)
                        data = next_line.split()
                        if len(data) < 5:
                            raise BadNetworkInput('Element type "fan" data point at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                        pts.append({'a1': self.handle_float(data[0]), 'a2': self.handle_float(data[1]),
                                    'a2': self.handle_float(data[3]), 'a4': self.handle_float(data[4]), 'mf': self.handle_float(data[4])})
                    obj['pts'] = pts
                    return obj
                # element name cpf upo prmin ftype
                elif data[2] == 'cpf':
                    if len(data) < 6:
                        raise BadNetworkInput('Element type "cpf" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    return {'input_type': InputType.ELEMENT, 'type': ElementType.CPF, 'name': data[1], 'upo': self.handle_float(data[3]),
                            'prmin': self.handle_float(data[4]), 'ftyp': self.handle_float(data[5])}
                # element name gfr a b
                elif data[2] == 'qfr':
                    if len(data) < 5:
                        raise BadNetworkInput('Element type "qfr" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    return {'input_type': InputType.ELEMENT, 'type': ElementType.QFR, 'name': data[1], 'a': self.handle_float(data[3]), 'b': self.handle_float(data[4])}
                # element name ckv dp0 coef
                elif data[2] == 'ckv':
                    if len(data) < 5:
                        raise BadNetworkInput('Element type "ckv" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
                    return {'input_type': InputType.ELEMENT, 'type': ElementType.CKV, 'name': data[1], 
                            'dp0': self.handle_float(data[3]), 'coeff': self.handle_float(data[4])}
                else:
                    raise BadNetworkInput('Element type "%s" not recognized' % data[2])
            elif data[0] == 'link':
                # link name node-l ht-l node-2 ht-2 element wind wpmod
                if len(data) < 8:
                    raise BadNetworkInput('Link at line %d has fewer than 8 fields' % self.line_number)
                if data[7] == 'null':
                    return {'input_type': InputType.LINK, 'name': data[1], 'node-1': data[2], 'ht-1': self.handle_float(data[3]), 
                            'node-2': data[4], 'ht-2': self.handle_float(data[5]), 'element': data[6]}
                if len(data) < 9:
                    raise BadNetworkInput('Link at line %d has fewer than 8 fields' % self.line_number)
                return {'input_type': InputType.LINK, 'name': data[1], 'node-1': data[2], 'ht-1': self.handle_float(data[3]), 
                        'node-2': data[4], 'ht-2': self.handle_float(data[5]), 'element': data[6],
                        'wind': data[7], 'wpmod': self.handle_float(data[8])}
            next_line = self.get_next()
        raise StopIteration