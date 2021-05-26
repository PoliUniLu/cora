

class ValuedVariable:
    
    def __init__(self, identificator, value):
        self._ident = identificator
        self._val = value
        
    def __hash__(self):
        return hash((self._ident, self._val))
    
    def __eq__(self, other):
        return self._ident == other._ident and self._val == other._val
    
    def __str__(self):
        return '{}({})'.format(self._ident, self._val)
    
    def __repr__(self):
        return str(self)

def bool_multiply(m):
    sets_to_multiply = []
    for row in m:
        s = {frozenset([ValuedVariable(i, v) 
                        ]) for i,v in enumerate(row)
                        if v != -1}
        sets_to_multiply.append(s)
    
    if len(sets_to_multiply) == 0:
        return set()
    
    res = sets_to_multiply[0]
    for k in range(1,len(sets_to_multiply)):
        tmp_res = {}
        for i,x in enumerate(res):
            for j,y in enumerate(sets_to_multiply[k]):
                new_element = x.union(y)
                if any(x.issubset(new_element) for x in tmp_res):
                    continue
                
                # Remove result elements which are subset of new element
                tmp_res = {x for x in tmp_res if not new_element.issubset(x)}
                tmp_res.add(new_element)
        res = tmp_res
    return res

def transform_to_raw_implicant(impl, levels):
    res = [frozenset(range(i)) for i in levels]
    for x in impl:
        res[x._ident] = frozenset([x._val])
    return tuple(res)


