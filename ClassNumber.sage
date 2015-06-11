import time

def qr(n,p):  #Legendre symbol
    if n == 0 or n % p == 0:
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
            
def kron_two(a,b): #kronecker symbol for p=2
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
        return kron_two(a,b)
        
    F = factor(b)
    F_lst = list(F)
    ks = kron_unit(a, F.unit())
    for p in F_lst:
        if p[0] == 2:
            ks = ks * (kron_two(a,b))^(p[1])
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

def E(x):
  series = 1
  for c in range(25, 0, -1):
    series = (2 * c + x - (c * (c + 1))/series)
    # print 2*c, pow(-1, c+1), c * (c + 1)
  series = 1 - 1/series
  series *= (e^-x)/x
  return series.numerical_approx()
  # return exp_integral_e1(x).numerical_approx()

def erfc(x):
  # series = 1
  # X = (x^2 - 0.5).numerical_approx()
  # for c in range(100, 0, -1):
  #   series = (2 * c + X - (c * (2 * c + 1)/2)/series)
  # series = (1 - 0.5/series).numerical_approx()
  # series *= e^(-x^2)/(x * sqrt(pi))
  # return series.numerical_approx()
  return (1 - erf(x)).numerical_approx()

def analyticClassNumberTheorem(d, units):
  for i in range(sqrt(d)+1, 1, -1):
    if d % (i * i) == 0:
      d //= (i * i)
  if d % 4 == 1:
    D = d
  else:
    D = 4 * d
  series = 0
  for r in range(1, int((D-1)/2)+1):
    series += (kron(D, r) * ln(sin((r * pi)/D))).numerical_approx()
  series /= -(ln(units[d])).numerical_approx()
  return series

def lam(D): #takes D = discriminant
    factor_lst = list(factor(D))
    for f in factor_lst:
        if f[0] % 4 == 3:
            return 2 #returns 2 if D has a prime divisor congruent to 3 (mod 4)
    return 1 #returns 1 if all prime divisors of D are congruent to 1 or 2 (mod 4)
    
def dist_prime(d): #takes d = fundamental discriminant
    return len(list(factor(d))) #returns the number of disctinct prime divisors of delta   

def L(D, R): #L-function to compute class number
  delta = D
  if D % 4 == 2 or D % 4 == 3:
      delta = 4 * delta
      
  t = dist_prime(delta) #number of distinct prime divisors
  s = lam(D) #takes value of 1 or 2
  l = ln((sqrt(pi/delta))/(R * 2^(t-s)))
  if l > 1:
      c = 6 - sqrt(27 - 2*l)
  elif l < 1:
      c = sqrt(15 + l) - 3  
  M = floor(c*sqrt(delta)/sqrt(pi)) + 1 #upper limit on sum to compute
  M *= 2 
  M = 2 * sqrt(D)
  C_m = 0
  i = 1
  while i <= M:
      C_m = C_m + kron(D,i) * (E(i^2 * pi/D) + erfc(i * sqrt(pi/D)) * sqrt(D)/i)
      i = i + 1
  return C_m

def slow_L(D): #L-function to compute class number
  delta = D
  if D % 4 == 2 or D % 4 == 3:
      delta = 4 * delta
      
  t = dist_prime(delta) #number of distinct prime divisors
  s = lam(D) #takes value of 1 or 2
  l = ln((sqrt(pi/delta))/(2^(t-s)))
  if l > 1:
      c = 6 - sqrt(27 - 2*l)
  elif l < 1:
      c = sqrt(15 + l) - 3  
  M = floor(c*sqrt(delta)/sqrt(pi)) + 1 #upper limit on sum to compute 
  M = 6 * sqrt(D)
  x,y,t = var('x, y, t')
  assume(x > 0)
  assume(y > 0)
  E = integral((e^(-t))/t, t, x, +Infinity)
  erfc = (2/sqrt(pi))*integral(e^(-t^2), t, y, +Infinity)
  C_m = 0
  i = 1
  while i <= M:
      C_m = C_m + kron(D,i) * (E(i^2 * pi/D) + erfc(i * sqrt(pi/D)) * sqrt(D)/i)
      i = i + 1
  return C_m

def ty_sum(d, units, stored):
  if d in stored:
    return stored[d]
  for i in range(sqrt(d)+1, 1, -1):
    if d % (i * i) == 0:
      d //= (i * i)
  if d % 4 == 1:
    D = d
  else:
    D = 4 * d
  return N((L(D, ln(units[d]))/(2 * ln(units[d])))).round()


stored = {}
# started 1:14
f = open("fundamental_units_(first_half).txt")
# f = open("medium_list")
# f = open("long_list")
# f = open("short_list")
# startTime = time.clock();
answers = open("faster_class_numbers(long_list_starting_at_300000)", "w")
unit_list = str.split(f.read(), "\n")
units = {}
for s in unit_list:
  components = str.split(s, " ", 1)
  # print components
  units[int(components[0])] = sage_eval(components[1])

print "made list"

step = 100
start = 300000
startTime = time.clock();
for j in range(1000):
  solns = ""
  for i in range(start + step*j, start + step*(j+1)):
    if sqrt(i) not in ZZ:
      # K = QuadraticField(i, 'x')
      calculated = int(ty_sum(i, units, stored))
      # if calculated != K.class_number():
      #   print str(i) + " wrong "
      #   solns += (str(i) + "\n")
      # else:
        # print str(i)
      solns += (str(i) + " " + str(calculated)) + "\n"
  answers.write(solns)
  answers.flush()  
  print(str(start + step*(j+1)) + " finished")
answers.close()
f.close()
print str(time.clock() - startTime) + " seconds"

# for d in units:
#   h = N((sqrt(d)/2)*(L(d)/ln(units[d]))).round()
#   print str(d) + " " + str(L(d).numerical_approx()) + " " + str(h)

# for d in units:
#   K = QuadraticField(d, 'x')
#   calculated = int(ty_sum(d, units, stored))
#   if calculated != K.class_number():
#     print str(d) + " wrong ", calculated, K.class_number()
#     # answers.write(str(d) + "\n")
#   else:
#     print str(d)
    # answers.write(str(d) + " " + str(calculated))

# answers.close()
# f.close()





# for e in fund_unit_lst:
#     h = N((sqrt(d)/2)*(L(d)/ln(e))).round()
#     h_lst.append(h)

# took 55.41 seconds
# f = open("fundamental_units(first_thousand).txt")
# startTime = time.clock();
# answers = open("class_numbers(first_thousand)", "w")
# unit_list = str.split(f.read(), "\n")
# units = {}
# for s in unit_list:
#   components = str.split(s, " ", 1)
#   # print components
#   units[int(components[0])] = sage_eval(components[1])

# print "made list"

# for i in units:
#   K = QuadraticField(i, 'x')
#   calculated = analyticClassNumberTheorem(i, units).round()
#   if calculated != K.class_number():
#     print str(i) + "wrong "
#   else:
#     print str(i)
#   answers.write(str(i) + " " + str(calculated) + "\n")

# f.close()
# print str(time.clock() - startTime) + " seconds"
# answers.close()

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


# took 9229 seconds
# step = 100
# start = 1
# f  =  open("fundamental_units(first_thousand)2.txt", "w")
# startTime = time.clock();
# for j in range(10):
#   solns = ""
#   for i in range(start + step*j, start + step*(j+1)):
#     if sqrt(i) not in ZZ:
#       answer = fundamentalUnit(i)
#       solns  += str(i) + " " + str(answer) + "\n"
#   f.write(solns)
#   f.flush()  
#   print(str(start + step*(j+1)) + " finished")
# f.close()
# print str(time.clock() - startTime) + " seconds"