import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import sympy as sp
l0 = sp.Symbol(r"\ell_{0}",Positive=True,Real=True)
I_r_m = sp.MatrixSymbol(r"{}^{\mathcal{I}}\mathbf{r}_m",3,1)
I_r_s = sp.MatrixSymbol(r"{}^{\mathcal{I}}\mathbf{r}_s", 3, 1)

I_S_matRot = sp.MatrixSymbol(r"{}_{\mathcal{S}}^{\mathcal{I}}\mathcal{R}", 3, 3)
I_M_matRot = sp.MatrixSymbol(r"{}_{\mathcal{M}}^{\mathcal{I}}\mathcal{R}", 3, 3)
S_b_k = sp.MatrixSymbol(r"{}^{\mathcal{S}}\mathbf{b}_k", 3, 1)
M_m_k = sp.MatrixSymbol(r"{}^{\mathcal{M}}\mathbf{m}_k", 3, 1)

f = I_r_m - I_r_s + I_M_matRot*M_m_k - I_S_matRot*S_b_k
sp.sqrt(f.T * f ) - l0