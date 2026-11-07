"""Exact identity checks used while deriving the saved formulas."""
import sympy as s
def simp(value): return s.factor(s.cancel(value))
class Identities:
    def __init__(self): self.names=[]
    def ok(self,name,value):
        if not bool(value): raise ArithmeticError(name)
        self.names.append(name)
    def zero(self,name,value):
        if isinstance(value,s.MatrixBase): values=list(value)
        elif isinstance(value,(list,tuple)): values=value
        else: values=[value]
        self.ok(name,all(simp(v)==0 for v in values))
