# Optional Colab/Jupyter cell source. Main reproducible program is odd_m_analysis.py.
from odd_m_analysis import group_data, diagonal_data
from test_odd_m import test_odd_m_cases
print('Exact cases verified:', test_odd_m_cases())
for family in ['M','N','P','O','W','Z']:
    for U,V in [(3,3),(3,2),(6,3),(6,2)]:
        print(family, U, V, group_data(family,U,V)['Gsize'])
