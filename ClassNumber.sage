def qr(n,p):  #Legendre symbol
    if n == 0:
        return 0
    else:
        if (n^((p-1)/2)) % p == 1:
            return 1
        else:
            return -1
            
def kron_unit(a,u): #kronecker symbol for units
    if u == 1:
      return 1
        
    elif u == -1:
        if a >= 0:
            return 1
            
        elif a < 0:
            return -1
            
def kron_two(a): #kronecker symbol for p=2
    if a % 2 == 0:
        return 0
            
    elif a % 8 == 1 or a % 8 == 7:
        return 1
            
    elif a % 8 == 3 or a % 8 == 5:
        return -1
        
def kron(a,b): #General kronecker symbol
    if b == 0:
        if a == 1 or a == -1:
            return 1
        else:
            return 0
    
    elif b == 1 or b == -1:
        return kron_unit(a,b)
            
    elif b == 2:
        return kron_two(b)
        
    F = factor(b)
    F_lst = list(F)
    ks = kron_unit(a, F.unit())
    for p in F_lst:
        if p[0] == 2:
            ks = ks * (kron_two(a))^(p[1])
        else:
            ks = ks * (qr(a,p[0]))^(p[1])
            
            
    return ks

# if item is in lst, returns its position. Otherwise returns -1
def in_list(lst, item):
  for x in range(len(lst)):
    # print(str(float(lst[x])) + " " + str(float(item)))
    if lst[x] == item:
      return x
  return -1

def gcd(a, b):
  if (a < 0):
    a *= -1
  if (b < 0):
    b *= -1

  while (b != 0):
    temp = b
    b = a % b
    a = temp

  return a


def reduce(n):
  # print(str(n[0]) + " " + str(n[1]) + " " + str(n[2]))
  divisor = gcd(n[0], n[1])
  divisor = gcd(divisor, n[2])
  # print(divisor)
  return (n[0]//divisor, n[1]//divisor, n[2]//divisor)

# returns the period and the continued fraction expansion of sqrt(n) in a tuple. If no repetition is found,
# returns -1 for the period
def rootContinuedFraction(n):
  alphas = []
  coeffs = []
  num = n
  root = sqrt(n)
  coeffs.append(int(root))
  alphas.append((0, 1, 1)) # store alphas in a triple (a, b, c) where alpha = (a + b*root(n))/c. This keeps the coefficients manageable
  n = (0, 1, 1)
  x = 0
  while(True):
    x+=1
    if n[1] == 0 and n[0]/n[2] in ZZ:
      coeffs.append(n[0]/n[2])
      return(-1, coeffs)
      break
    else:
      # alpha_{n+1} = 1/(a_n + alpha_n). The expression below is what happens when you multiply above and below by the conjugate of a_n + alpha_n
      n = (n[2]*(coeffs[x-1]*n[2] - n[0]), n[1]*n[2], n[1]**2*num - (coeffs[x-1]*n[2]-n[0])**2)
      
      # reduce the fraction by dividing out by the gcd of the numerator and denominator
      n = reduce(n)
      
      # check if we have computed this alpha already. If so, we will repeat from this point on. The distance between then and now is the period
      pos = in_list(alphas, n)
      
      if pos != -1:
        return (len(coeffs)-pos, coeffs)

      # otherwise, keep going
      alphas.append(n)

      # a is the floor of
      a = int((n[0] + n[1] * root)/n[2])
      coeffs.append(int(a))

  return (-1, coeffs)

# returns the kth convergent of fraction frac. Frac is a tuple. The first component is the period of the fraction, the second is the continued fraction upon to the first repetiti
def fromContinuedFraction(frac, k):
  if k == 0:
    return (frac[1][0], 1)

  # record when the continued fraction begins to repeat
  startOfCycle = len(frac[1])-frac[0]

  # now, we will store k places of the continued fraction
  explicitFrac = []
  for x in range(k+1):
    # until the record of the fraction ends, store the numbers normally
    if (x < len(frac[1])):
      explicitFrac.append(frac[1][x])

    # once the fraction begins to repeat, we have to look to the end of frac to find the right numbers
    else:
      explicitFrac.append(frac[1][((x - len(frac[1]))%frac[0] + startOfCycle)])

  # store the numerator and denominator separately
  numerator = explicitFrac[k]
  denominator = 1
  for x in range(k-1, -1, -1):
    old_numerator = numerator

    # we perform fraction addition a + 1/(b/c) = (ba + c)/b
    numerator = explicitFrac[x] * numerator + denominator
    denominator = old_numerator

    # reduce the fraction by dividing by the gcd
    divisor = gcd(numerator, denominator)
    numerator /= divisor
    denominator /= divisor
  return (numerator, denominator)


# finds the fundamental solution of the positive Pell equation x^2 - dy^2 = 1 using continued fractions. Returns this solution as a tuple (x, y)
def fundamentalSolution(d):
  # find the continued fracgtion expansion of d. This returns a tuple. The first component is the period of the continued fraction. The second component is the fraction itself
  continued_frac = rootContinuedFraction(d)
  period = continued_frac[0]

  return fromContinuedFraction(continued_frac, period - 1)

  # if period % 2 == 0:
  #   return fromContinuedFraction(continued_frac, period - 1)
  # else :
  #   return fromContinuedFraction(continued_frac, 2 * period - 1)

# finds the fundamental unit of the quadratic field Q(sqrt(d))
def fundamentalUnit(d):
  if d % 8 != 5:
    soln = fundamentalSolution(d)
    return (soln[0] + sqrt(d) * soln[1])
  else:
    soln=fundamentalSolution(d)
    # print(str(soln[0]) + " " + str(soln[1]))
    for a in range(1, 2 * pow(soln[0], 1/3) + 1):
      b = sqrt((8*soln[0] - pow(a, 3))/(3*a*d))
      if b in ZZ and 3 * pow(a, 2) * b + pow(b, 3) * d == 8 * soln[1]:
          return (a/2 + b/2 * sqrt(d))
    return (soln[0] + soln[1] * sqrt(d))

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

def analyticClassNumberTheorem(d, units):
  # print "original d: " + str(d)
  for i in range(sqrt(d)+1, 1, -1):
    # print str(i) + " " +  str(d % i * i)
    if d % (i * i) == 0:
      d //= (i * i)
  # print "d: " + str(d)
  if d % 4 == 1:
    D = d
  else:
    D = 4 * d
  series = 0
  for r in range(1, int((D-1)/2)+1):
    series += (kronecker_symbol(D, r) * ln(sin((r * pi)/D))).numerical_approx()
  series /= -(ln(units[d])).numerical_approx()
  return series


#started 11:54
# f = open("fundamental_units(first_thousand).txt")
# f = open("short_list")
# f = open("longer_list")
f = open("fundamental_units_cohen_all.txt")
answers = open("class_numbers(all)", "w")
unit_list = str.split(f.read(), "\n")
units = {}
for s in unit_list:
  components = str.split(s, " ", 1)
  print components
  units[int(components[0])] = sage_eval(components[1])

print "made list"

for i in units:
  K = QuadraticField(i, 'x')
  calculated = analyticClassNumberTheorem(i, units).round()
  if calculated != K.class_number():
    print str(i) + "wrong "
  else:
    print str(i)
  answers.write(str(i) + " " + str(calculated) + "\n")

f.close()
answers.close()
  # print str(i) + " | correct: " + str(K.class_number()) + " | calculated: " + str((analyticClassNumberTheorem(i, units)).round())
    

# started 9:55
# finished 10:50
# step = 1000
# start = 3
# f  =  open("fundamental_units_cohen_all.txt", "w")
# for j in range(1000):
#   solns = ""
#   for i in range(start + step*j, start + step*(j+1)):
#     if sqrt(i) not in ZZ:
#       answer = cohenFundamentalUnitSlower(i)
#       solns  += str(i) + " " + str(answer) + "\n"
#   f.write(solns)
#   f.flush()  
#   print(str(start + step*(j+1)) + " finished")
# f.close()


# for i in range(15):



# for x in range(25):
#   for y in range(25):
#     if kron(x, y) != kronecker_symbol(x, y):
#       print("(" + str(x) + " , " + str(y) + ") | kron: " + str(kron(x, y)) + " | kronecker: " + str(kronecker_symbol(x, y)))
# print "done:"