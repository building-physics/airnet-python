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

#class ElementType(enum.Enum):
#    PLR = 0
#    DWC = 1
#    DOR = 2
#    CFR = 3
#    FAN = 4
#    CPF = 5
#    QFR = 6
#    CKV = 7

class Reader:
    def __init__(self, fp, floats=True):
        self.fp = fp
        self.title = None
        self.line_number = 0
        self.handle_float = self.handle_float_as_string
        if floats:
            self.handle_float = self.handle_float_as_float

    def get_next(self):
        next_line = next(self.fp)
        self.line_number += 1
        while(next_line.strip() == ''):
            next_line = next(self.fp)
            self.line_number += 1
        return next_line
    
    def __iter__(self):
        return self
    
    def handle_float_as_string(self, type, field, value):
        try:
            float(value)
        except ValueError:
            raise BadNetworkInput('Value "%s" for field "%s" of "%s" at line %d in not a float'
                                  % (value, field, type, self.line_number))
        return value
    
    def handle_float_as_float(self, type, field, value):
        try:
            return float(value)
        except ValueError:
            raise BadNetworkInput('Float conversion of "%s" failed in "%s" for field "%s" at line %d'
                                  % (value, type, field, self.line_number))

    def __next__(self):
        next_line = self.get_next()
        while not next_line.startswith('*'):
            next_line = next_line.lstrip()
            item_type, sep, rest = next_line.partition(' ')
            if item_type == 'title':
                if self.title:
                    raise BadNetworkInput('Found additional title at line %d' % self.line_number)
                else:
                    self.title = rest.strip()
                    return {'input_type': InputType.TITLE, 'title': self.title}
            
            if item_type == 'node':
                name, sep, input_string = rest.lstrip().partition(' ')
                return self.read_node(name, input_string.lstrip())

            elif item_type == 'element':
                name, sep, input_string = rest.lstrip().partition(' ')
                element_type, sep, input_string = input_string.lstrip().partition(' ')
                method_name = "read_%s" % element_type
                try:
                    result = getattr(self, method_name)(name, input_string.lstrip())
                except AttributeError:
                    raise BadNetworkInput('Element type "%s" not recognized' % element_type)
                return result
            
            elif item_type == 'link':
                name, sep, input_string = rest.lstrip().partition(' ')
                return self.read_link(name, input_string.lstrip())
            next_line = self.get_next()
        raise StopIteration
    
    def read_node(self, name, input_string):
        data = input_string.split()
        # node name type ht temp pres
        if len(data) < 3:
            raise BadNetworkInput('Node at line %d has fewer than 5 fields' % self.line_number)
        # node name type ht temp pres
        if data[0] not in ['v', 'c', 'a']:
            raise BadNetworkInput('Node type "%s" at line %d is unrecognized, must be "v", "c", or "a"' % (data[0], self.line_number))
        if data[0] == 'v':
            return {'input_type': InputType.NODE, 'name': name, 'type': data[0],
                    'ht': self.handle_float('node', 'ht', data[1]), 
                    'temp': self.handle_float('node', 'temp', data[2])}
        else:
            if len(data) < 4:
                raise BadNetworkInput('Node at line %d has fewer than 6 fields' % self.line_number)
            return {'input_type': InputType.NODE, 'name': name,
                    'type': data[0], 'ht': self.handle_float('node', 'ht', data[1]), 
                    'temp': self.handle_float('node', 'temp', data[2]),
                    'pres': self.handle_float('node', 'pres', data[3])}
        

    def read_link(self, name, input_string):
        data = input_string.split()
        # link name node-l ht-l node-2 ht-2 element wind wpmod
        if len(data) < 6:
            raise BadNetworkInput('Link at line %d has fewer than 8 fields' % self.line_number)
        if data[5] == 'null':
            return {'input_type': InputType.LINK, 'name': name, 'node-1': data[0],
                    'ht-1': self.handle_float('link', 'ht-l', data[1]), 
                    'node-2': data[2], 'ht-2': self.handle_float('link', 'ht-2', data[3]), 'element': data[4]}
        if len(data) < 7:
            raise BadNetworkInput('Link at line %d has fewer than 9 fields' % self.line_number)
        return {'input_type': InputType.LINK, 'name': name, 'node-1': data[0], 
                'ht-1': self.handle_float('link', 'ht-l', data[1]), 
                'node-2': data[2], 'ht-2': self.handle_float('link', 'ht-2', data[3]), 'element': data[4],
                'wind': data[5], 'wpmod': self.handle_float('link', 'wpmod', data[6])}
    
    def read_plr(self, name, input_string):
        # element name plr init lam turb expt
        data = input_string.split()
        if len(data) < 4:
            raise BadNetworkInput('Element type "plr" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)+3))
        return {'input_type': InputType.ELEMENT, 'type': 'plr', 'name': name,
                'init': self.handle_float('plr', 'init', data[0]), 'lam': self.handle_float('plr', 'lam', data[1]),
                'turb': self.handle_float('plr', 'turb', data[2]), 'expt': self.handle_float('plr', 'expt', data[3])}

    def read_dwc(self, name, input_string):
        # element name dwc len dh area rgh
        #         tdlc lflc ldlc init
        data = input_string.split()
        if len(data) < 4:
            raise BadNetworkInput('Element type "dwc" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)+3))
        obj = {'input_type': InputType.ELEMENT, 'type': 'dwc', 'name': name,
               'len': self.handle_float('dwc', 'len', data[0]), 
               'dh': self.handle_float('dwc', 'dh', data[1]), 'area': self.handle_float('dwc', 'area', data[2]),
               'rgh': self.handle_float('dwc', 'rgh', data[3])}
        try:
            next_line = self.get_next()
        except StopIteration:
            raise BadNetworkInput('Element type "dwc" at line %d has only one line and cannot be a legal element' % self.line_number)
        data = next_line.split()
        if len(data) < 4:
            raise BadNetworkInput('Element type "dwc" with second line at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
        obj['tdlc'] = self.handle_float('dwc', 'tdlc', data[0]) 
        obj['lflc'] = self.handle_float('dwc', 'lfdc', data[1])
        obj['ldlc'] = self.handle_float('dwc', 'ldlc', data[2]) 
        obj['init'] = self.handle_float('dwc', 'init', data[3])
        return obj
    def read_dor(self, name, input_string):
        # element name dor init lam turb expt
        #         dtmin ht wd cd
        data = input_string.split()
        if len(data) < 4:
            raise BadNetworkInput('Element type "dor" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)+3))
        obj = {'input_type': InputType.ELEMENT, 'type': 'dor', 'name': name,
               'init': self.handle_float('dor', 'init', data[0]), 'lam': self.handle_float('dor', 'lam', data[1]),
               'turb': self.handle_float('dor', 'turb', data[2]), 'expt': self.handle_float('dor', 'expt', data[3])}
        try:
            next_line = self.get_next()
        except StopIteration:
            raise BadNetworkInput('Element type "dor" at line %d has only one line and cannot be a legal element' % self.line_number)
        data = next_line.split()
        if len(data) < 4:
            raise BadNetworkInput('Element type "dor" with second line at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
        obj['dtmin'] = self.handle_float('dor', 'dtmin', data[0]) 
        obj['ht'] = self.handle_float('dor', 'ht', data[1])
        obj['wd'] = self.handle_float('dor', 'wd', data[2]) 
        obj['cd'] = self.handle_float('dor', 'cd', data[3])
        return obj
    
    def read_cfr(self, name, input_string):
        # element name cfr flow
        string = input_string.partition(' ')[0].strip()
        if not string:
            raise BadNetworkInput('Element type "cfr" at line %d has only 3 fields and cannot be a legal element' % self.line_number)
        return {'input_type': InputType.ELEMENT, 'type': 'cfr', 'name': name,
                'flow': self.handle_float('cfd', 'flow', string)}
    
    def read_fan(self, name, input_string):
        # element name fan init lam turb expt
        #         rdens fdf sop off nr mfl
        #         all al2 al3 a14 mf2
        #         a2l a22 a23 a24 mf3
        #         ... ... ... ... mfn
        data = input_string.split()
        if len(data) < 7:
            raise BadNetworkInput('Element type "fan" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)+3))
        obj = {'input_type': InputType.ELEMENT, 'type': 'fan', 'name': name,
               'init': self.handle_float('fan', 'init', data[0]), 'lam': self.handle_float('fan', 'lam', data[1]),
               'turb': self.handle_float('fan', 'turb', data[2]), 'expt': self.handle_float('fan', 'expt', data[3])}
        try:
            next_line = self.get_next()
        except StopIteration:
            raise BadNetworkInput('Element type "fan" at line %d has only one line and cannot be a legal element' % self.line_number)
        data = next_line.split()
        if len(data) < 6:
            raise BadNetworkInput('Element type "dwc" with second line at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
        obj['rdens'] = self.handle_float('fan', 'rdens', data[0]) 
        obj['fdf'] = self.handle_float('fan', 'fdf', data[1])
        obj['sop'] = self.handle_float('fan', 'sop', data[2]) 
        obj['off'] = self.handle_float('fan', 'off', data[3])
        nr = int(data[4]) 
        obj['mf1'] = self.handle_float('fan', 'mfl', data[5])
        pts = []
        for i in range(nr):
            try:
                next_line = self.get_next()
            except StopIteration:
                raise BadNetworkInput('Element type "fan" at line %d has insufficient lines of data and cannot be a legal element' % self.line_number)
            data = next_line.split()
            if len(data) < 5:
                raise BadNetworkInput('Element type "fan" data point at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)))
            pts.append({'a1': self.handle_float('fan', 'a1%d' % (i+1), data[0]),
                        'a2': self.handle_float('fan', 'a2%d' % (i+1), data[1]),
                        'a3': self.handle_float('fan', 'a3%d' % (i+1), data[3]),
                        'a4': self.handle_float('fan', 'a4%d' % (i+1), data[4]), 
                        'mf': self.handle_float('fan', 'mf%d' % (i+1), data[4])})
        obj['pts'] = pts
        return obj
    
    def read_cpf(self, name, input_string):
        # element name cpf upo prmin ftype
        data = input_string.split()
        if len(data) < 3:
            raise BadNetworkInput('Element type "cpf" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)+3))
        return {'input_type': InputType.ELEMENT, 'type': 'cpf', 'name': name,
                'upo': self.handle_float('cpf', 'upo', data[0]),
                'prmin': self.handle_float('cpf', 'prmin', data[1]), 'ftyp': self.handle_float('cpf', 'ftyp', data[2])}
    
    def read_qfr(self, name, input_string):
        # element name qfr a b
        data = input_string.split()
        if len(data) < 2:
            raise BadNetworkInput('Element type "qfr" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)+3))
        return {'input_type': InputType.ELEMENT, 'type': 'qfr', 'name': name, 
                'a': self.handle_float('qfr', 'a', data[0]), 'b': self.handle_float('qfr', 'b', data[1])}
    
    def read_ckv(self, name, input_string):
        # element name ckv dp0 coef
        data = input_string.split()
        if len(data) < 2:
            raise BadNetworkInput('Element type "ckv" at line %d has only %d fields and cannot be a legal element' % (self.line_number, len(data)+3))
        return {'input_type': InputType.ELEMENT, 'type': 'ckv', 'name': name, 
                'dp0': self.handle_float('ckv', 'dp0', data[0]), 'coef': self.handle_float('ckv', 'coef', data[1])}
