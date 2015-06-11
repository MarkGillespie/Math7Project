import time

#page 270
def cohenFundamentalUnitSlower(D):
  d = int(sqrt(D))

  if (d % 2 == D % 2):
    b = d
  else:
    b = d-1

  # uses continued fractions to expand (b + sqrt D)/2

  # alpha_i, the real number used in the continued fraction algorithm to find the ith term a_i = int(alpha_i),
  # is here given by alpha_i = (p_i + d)/q_i
  # the i subscript means this is true after iterating i times.
  # note that alpha_0 = (b + sqrt D)/2, the fraction we are expanding
  p = b
  q = 2

  # we look at the sequence a_i, b_i where  a_i/b_i is the ith convergent

  # a_{n-1} and a_n
  u1 = -b
  u2 = 2

  # b_{n-1} and b_n
  v1 = 1
  v2 = 0

  loop = True

  while loop:
    # A is the ith term in the continued fraction expansion of (b + sqrt D)/2
    A = int((p + d)/q)

    # update p and q to be p_{i+1} and q_{i+1}
    p = A * q - p
    q = (D - pow(p, 2))/q


    # u, v updated using the matrix rule
    # 
    #   u2' v2'   A  1     u2  v2
    #   u1' v1' = 1  0  x  u1  v1
    #  
    # this is a compact way of writing the recursion used to calculate convergents
    t = A * u2 + u1
    u1 = u2
    u2 = t

    t = A * v2 + v1
    v1 = v2
    v2 = t

    # stop looping when q = 2 and p = b (mod 2)
    loop = not (q == 2 and p % 2 == b % 2)

  # cohen proves that when q = 2 and p = b (mod 2), this convergent is a solution
  u = abs(u2)
  v = abs(v2)
  return(u/2 +  v/2 * sqrt(D))

units = {}

step = 1000
start = 3
f  =  open("fundamental_units_cohen_square.txt", "w")

# took 3204 seconds
startTime = time.clock()
for j in range(1000):
  solns = ""
  for d in range(start + step*j, start + step*(j+1)):
    if sqrt(d) not in ZZ:
      disc = d
      for i in range(sqrt(disc)+1, 1, -1):
        if disc % (i * i) == 0:
          disc //= (i * i)
      # print disc, d
      if disc == d:
        answer = cohenFundamentalUnitSlower(d)
        units[d] = answer
        solns  += str(d) + " " + str(answer) + "\n"
  f.write(solns)
  f.flush()  
  print(str(start + step*(j+1)) + " finished")
f.close()
print time.clock()-startTime