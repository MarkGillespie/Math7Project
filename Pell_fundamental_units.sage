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

# took 9229 seconds
step = 100
start = 1
f  =  open("fundamental_units(first_thousand)2.txt", "w")
startTime = time.clock();
for j in range(10):
  solns = ""
  for i in range(start + step*j, start + step*(j+1)):
    if sqrt(i) not in ZZ:
      answer = fundamentalUnit(i)
      solns  += str(i) + " " + str(answer) + "\n"
  f.write(solns)
  f.flush()  
  print(str(start + step*(j+1)) + " finished")
f.close()
print str(time.clock() - startTime) + " seconds"