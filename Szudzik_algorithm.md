# Site Mapping with Szudzik Algorithm 

```import numpy as np

def szudzik(a, b):
    if(a >= b):
        return a * a + a + b
    else:
        return a + b * b

def combine(i, j, k):
    return szudzik(szudzik(i, j), k)


def szudzik_mapping(i_max, j_max, k_max):
    mapping = {}
    
    for i in np.arange(0, i_max + 1): 
        for j in np.arange(0, j_max + 1):
            for k in np.arange(0, k_max + 1):
                mapping[combine(i, j, k)] = (i, j, k)
    
    return mapping
```