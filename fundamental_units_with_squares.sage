f = open("fundamental_units_cohen_square.txt")
answers = open("fundamental_units_cohen_all(fixed)", "w")
unit_list = str.split(f.read(), "\n")
units = {}

for s in unit_list:
  components = str.split(s, " ", 1)
  units[int(components[0])] = sage_eval(components[1])

print "made list"

for i in range(1000):
  stream = ""
  for j in range(1000):
    k = 1000 * i + j

    if sqrt(k) not in ZZ:
      d = k

      # d is the squarefree part of k
      for f in range(sqrt(d)+1, 1, -1):
        if d % (f * f) == 0:
          d //= (f * f)

      stream += str(k) + " " + str(units[d]) + "\n"
  answers.write(stream)
  answers.flush()
  print "done with ", 1000 * (i + 1)

answers.close()