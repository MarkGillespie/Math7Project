import time

#page 270
def cohenFundamentalUnitSlower(D):
  d = int(sqrt(D))
  if (d % 2 == D % 2):
    b = d
  else:
    b = d-1

  u1 = -b
  u2 = 2
  v1 = 1
  v2 = 0
  p = b
  q = 2
  loop = True

  while loop:
    # print (str(D) + " " + str(p) + " " + str(q))
    A = int((p + d)/q)
    p = A * q - p
    q = (D - pow(p, 2))/q
    t = A * u2 + u1
    u1 = u2
    u2 = t
    t = A * v2 + v1
    v1 = v2
    v2 = t

    loop = not (q == 2 and p % 2 == b % 2)

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