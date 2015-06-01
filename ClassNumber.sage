# Jacobi Symbol (n/p)
def qr(n, p):
    # base cases
    if n == 0 or n == p or mod(n, p) == 0:
        return 0
    # positive base case: 1 is a square
    elif n == 1:
        return 1
    # special case for n == 2
    # we proved in class that (2, p) = 1 iff p is congruent to plus or minus 1 modulo 8
    elif n == 2:
        if mod(p, 8) == 1 or mod(p, 8) == 7:
            return 1
        else:
            return -1
    # special case for p = 2
    # 1 is the only quadratic residue mod 2, so check if n is congruent to 1 mod 2
    elif p == 2:
        return 1
    # recursive case
    else:
        #only use quadratic reciprocity if p,q odd. Otherwise split by dividing one by 2
        if mod(p, 2) == 1:
            if mod(n, 2) == 1:
                if mod(n, 4) == 1 or mod(p, 4) == 1:
                    return qr(p % n,  n)
                else:
                    return -1*qr(p % n,  n)
            else:
                return qr(2, p) * qr(n/2, p)
        else:
            return qr(n, 2) * qr(n, p/2)

# if item is in lst, returns its position. Otherwise returns -1
def in_list(lst, item):
  for x in range(len(lst)):
    # print(str(float(lst[x])) + " " + str(float(item)))
    if lst[x] == item:
      return x
  return -1

def gcd(a, b):

  if (a < 0):
    return gcd(-a, b)
  if (b < 0):
    return gcd(a, -b)
  if b == 0:
    return a
  return gcd(b, a%b)

def reduce(n):
  # print(str(n[0]) + " " + str(n[1]) + " " + str(n[2]))
  divisor = gcd(n[0], n[1])
  divisor = gcd(divisor, n[2])
  # print(divisor)
  return (n[0]//divisor, n[1]//divisor, n[2]//divisor)

# returns the period and the continued fraction in a tuple. If no repetition is found,
# returns -1 for the period
def rootContinuedFraction(n, k):
  alphas = []
  coeffs = []
  num = n
  root = n**(0.5)
  coeffs.append(int(root))
  alphas.append((0, 1, 1)) # store alphas in a triple (a, b, c) where alpha = (a + b*root(n))/c. This keeps the coefficients manageable
  print(root)
  n = (0, 1, 1)
  for x in range(1,k):
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

print(rootContinuedFraction(31, 50))
# print(gcd(1256, 728))