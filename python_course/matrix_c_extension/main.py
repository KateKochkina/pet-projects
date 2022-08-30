import matrix

m = matrix.Matrix([[1, 2], [3, 4], [5, 6]])
n = matrix.Matrix([[1, 2, 3], [4, 5, 6]])

print(m)
print(m.sum(m))
print(m.mul_scalar(5))
print(m.div(2))
print(m.mul(n))
print(m.transp())
print(5 in m)
print(7 in m)
print(m[(0,1)])
