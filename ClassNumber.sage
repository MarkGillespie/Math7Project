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
  x,y,t = var('x, y, t')
  assume(x > 0)
  assume(y > 0)
  E = integral((e^(-t))/t, t, x, +Infinity)
  erfc = (2/sqrt(pi))*integral(e^(-t^2), t, y, +Infinity)
  while i <= M:
      term = kron(D,i) * (E(i^2 * pi/D) + erfc(i * sqrt(pi/D)) * sqrt(D)/i)
      C_m = C_m + term
      i = i + 1
  print D, C_m.numerical_approx()
  return C_m

def ty_sum(d, units, stored):
  q = d
  if d in stored:
    return stored[d]
  for i in range(sqrt(d)+1, 1, -1):
    if d % (i * i) == 0:
      d //= (i * i)
  if d % 4 == 1:
    D = d
  else:
    D = 4 * d
  s = N(L(D, ln(units[d])))
  print q, s , N(2 * ln(units[d]))
  return (s/(2 * ln(units[d]))).round()



stored = {}
# started 1:14
# f = open("fundamental_units_(first_half).txt")
# f = open("medium_list")
f = open("long_list")
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
start = 1#300000
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
#       solns += (str(i) + " " + str(calculated)) + "\n"
#   answers.write(solns)
#   answers.flush()  
#   print(str(start + step*(j+1)) + " finished")
# answers.close()
# f.close()
# print str(time.clock() - startTime) + " seconds"

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