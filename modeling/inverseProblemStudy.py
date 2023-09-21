import sympy as sp
from  common_functions.QuatSymbolic import QuaternionSymbolic

q = QuaternionSymbolic()

def screw(v):
    q_x = sp.Matrix(sp.zeros(3, 3))
    q_x[0, 1] = -v[-1]
    q_x[1, 0] = v[-1]
    q_x[0, 2] = v[-2]
    q_x[2, 0] = -v[-2]
    q_x[1, 2] = -v[-3]
    q_x[2, 1] = v[-3]
    return q_x

m = sp.Matrix([[sp.symbols("m_1")],
               [sp.symbols("m_2")],
               [sp.symbols("m_3")]])

q.Q().T @ q.Q() @ q.quatRot(normalized=1)


df_dq = sp.Matrix([q.quat[0]* m.T - m.T * q.screw(),
                   (q.vec().T*m)[0,0]*sp.eye(3)+m*q.vec().T-q.vec()*m.T+q.quat[0]*screw(m)])

q.quatRot(normalized=True)

t = q.Q().T @ df_dq@q.quatRot(normalized=True).T


I = sp.diag(sp.symbols(r'I_x'), sp.symbols(r'I_x'), sp.symbols(r'I_x'))**-1

# t = df_dq * q.quatRot(normalized=1).T
for lin in range(3):
    for col in range(3):
        t[lin,col] = t[lin,col].expand().simplify()
        t[lin,col] = t[lin, col].subs(q.quat[0]**2, 1 - (q.vec().T*q.vec())[0])
        t[lin, col] = t[lin, col].simplify()
        t[lin, col] = t[lin, col].factor(deep=True)
t
t[2,0].expand()


print(sp.latex(t))

q.screw().T @ q.screw()


t,f = sp.symbols(r't, f')
A = sp.symbols(r'A',positive=True)

p = A*(1-sp.exp(-t-A))
v = sp.diff(p,t)
a = sp.diff(v,t)
a
p
v
a

a.subs(t,0)

a = A*(1-sp.exp(-A*t)).collect(A)
v = sp.integrate(a,t).collect(A)
p = sp.integrate(v,t).collect(A)
a
v
p


p.subs(t, 0)

Rb = sp.MatrixSymbol(r'{}^{I}_{B}R',3,3)
Rm = sp.MatrixSymbol(r'{}^{I}_{m}R', 3, 3)
x = sp.MatrixSymbol(r'\chi',3,1)


b = sp.MatrixSymbol('b',3,1)
m = sp.MatrixSymbol('m',3,1)



f = x + (Rm * m ) - (Rb * b)
ff = (f.T * f).expand()






with sp.assuming(sp.Q.orthogonal(Rb)) , sp.assuming(sp.Q.orthogonal(Rm)):
    ff2 = sp.refine(ff)
    ff2


