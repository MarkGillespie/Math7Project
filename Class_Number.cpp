#include<Class_Number.h>

//Jacobi Symbl (n/p)
int qr(int n, int p){
  //positive base case: 1 is a square
  if (n == 1 || p == 1) {
    return 1;
  //base cases
  } else if (n == 0 || n == p || n % p == 0 || p % n == 0) {
    return 0;
  //special case for n == 2
  //we proved in class that (2, p) = 1 iff p is congruent to plus or minus 1 modulo 8
  } else if (n == 2) {
    if (p % 8 == 1 || p % 8 == 7) {
      return 1;
    } else {
      return -1;
    }
  //special case for p = 2
  //1 is the only quadratic residue mod 2, so check if n is congruent to 1 mod 2
  } else if (p == 2) {
    if (n % 8 == 1 || n % 8 == 7) {
      return 1;
    } else {
      return -1;
    }
  //recursive case
  } else {
    // only use quadratic reciprocity if p,q odd. Otherwise split by dividing one by 2
    if (p % 2 == 1) {
      if (n % 2 == 1) {
        if (n % 4 == 1 || p % 4 == 1) {
          return qr(p % n,  n);
        } else {
          return -1*qr(p % n,  n);
        }
      } else {
        return qr(2, p) * qr(n/2, p);
      }
    } else {
      return qr(n, 2) * qr(n, p/2);
    }
  }
}

double gamma_inc(double a, double z){
  double s = 1;
  for (int i = 100; i > 0; i--){
    s = 2 * i - 1 - a + z + (double)(i * (a - i))/s;
  }
  s = pow(z, a) * pow(e,-z)/s;
  return s;
}

double erfc(double x) {
  double s = 1;
  double X = pow(x, 2) - 0.5;

  for (int i = 100; i > 0; i--) {
    s = (double)2 * i + X - (i * (2 * i + 1)/(double) 2)/s;
  }
  s = 1 - 0.5/s;
  s *= pow(e, -pow(x, 2))/(x * sqrt(pi));
  return s;
}

double E(double x) {
  double series = 1;
  for (int c = 100; c > 0; c--) {
    series = (2 * c + x - (c * (c + 1))/series);
  }
  series = 1 - 1/series;
  series *= pow(e,-x)/x;
  return series;
}

// int get_class_number(int d, double R) {

//   double logd = log(d);
//   double p1 = sqrt(d * logd) / (2 * pi);
//   double p2 = 1 - 2 * R / logd;
//   if (pow(p2, 2) > 2 / logd)
//     p1 *= p2;

//   int n = round(p1);
//   double p4 = pi/d;

//   double rootd = sqrt(d);
//   double p5 = (double)1 - 1 / sqrt(pi) * gamma_inc(0.5, p4);
//   double S = rootd * p5 + E(p4);
//   int k;

//   for (int i = 2; i <= n; i++) {
//     k = qr(d, i);
//     if (k == 0)
//       continue;
//     p2 = p4 * sqrt(i);
//     p5 = 1 - 1/sqrt(pi) * gamma_inc(0.5, p2);
//     if (k > 0) {
//       S += (p2 + rootd * p5 / (double)i);
//     } else {
//       S -= (p2 + rootd * p5 / (double)i);
//     } 
//   }
//   std::cout<<d<<" "<<S<<" "<<2*R<<std::endl;
//   S /= (double)(2 * R);
//   return round(S);
// }

double L(int D, double R) {

  // double logd = log(D);
  // double p1 = sqrt(D * logd) / (2 * pi);
  // double p2 = 1 - 2 * R / logd;
  // if (pow(p2, 2) > 2 / logd)
  //   p1 *= p2;

  // int n = round(p1) + 5;
  int n = 4 * (int)sqrt(D);
  double s = 0;
  int k;
  for (int i = 1; i <= n; i++) {
    k = qr(D, i);
    if (k != 0) {
      double term = (double)k * (E(pow(i, 2) * pi/((double) D)) + erfc(i * sqrt(pi/( (double) D))) * sqrt(D)/((double) i));
      s += term;
      // std::cout<<"i: "<<i<<" s: "<<s<<" term: "<<term<<std::endl;
    }
  }
  // std::cout<<"D: "<<D<<" "<<s<<std::endl;
  return s;
} 

int get_class_number(int D, double R) {
  // std::cout<<"R: "<<R;
  return round(L(D, R)/(2 * R));
}

bool isSquare(int n) {
  int root = sqrt(n);
  return n == root * root;
}

void fill_units(std::string filename, int length, double(arr)[]) {
  std::ifstream myFile;
  myFile.open(filename);
  std::string line;

  if (myFile.is_open()) {
    int counter = 2;
    while (getline(myFile, line) && counter < length) {
      int split_pos = line.find(" ");
      int index = std::stoi(line.substr(0, split_pos));
      double unit = std::stod(line.substr(split_pos+1));
      // std::cout<<line.substr(split_pos+1)<<" "<<unit<<std::endl;
      arr[index] = unit;
      counter++;
      if (isSquare(counter)) {
        counter++;
      }
    }
  }
  myFile.close();
}



int main(int argc, char **argv) {
  // std::cout<< class_number(3, log(sqrt(3) + 2))<<std::endl;
  
  std::ofstream output ("c_class_numbers.txt");
  if (!output.is_open()) {
    std::cout<<"closed"<<std::endl;
  }

  int input_len = 1000000;
  double units[1000000] = {0};
  int class_nums[1000000] = {0};
  // std::cout<<units[0]<<std::endl;
  fill_units("short_list_doubles", input_len, units);

  int step = 2;
  int giant_step = input_len/step;
  for (int i = 0; i < giant_step; i++) {
    std::string stream = "";
    for (int j = 0; j < step; j++) {
      int position = i * step + j;
      if (!isSquare(position)) {
        int d = position;
        int class_number = 0;
        for (int divisor = (int)(sqrt(d)); divisor > 1; divisor--) {
          if (d % (int)pow(divisor, 2) == 0) {
            d /= divisor;
            d /= divisor;
          }
        }
        if (d != position) {
          class_number = class_nums[d];
        } else {
          int D = d;
          if (d % 4 != 1) {
            D *= 4;
          }
          class_number = get_class_number(D, log(units[d]));
        }
        class_nums[i] = class_number;
        stream += std::to_string(position) + " " + std::to_string(class_number) + (std::string)"\n";
      }      
    }
    // std::cout<<stream<<std::endl;
    output<<stream;
    std::cout<<"finished with "<<step * (i + 1)<<std::endl;
    
  }
  output.close();

  
  // for (int i = 1; i < 10; i++) {
  //   for (int j = 1; j < 10; j++) {
  //     std::cout<<qr(i, j)<<" ";
  //   }
  //   std::cout<<std::endl;
  // }
  // std::cout<<qr(7, 3);
  
  // std::cout<<qr(7, 1)<<std::endl;
  // for (int i = 2; i < 25; i++) {
  //   // std::cout<<"gamma: " <<gamma_inc(0.5, i)<<" E: "<<E(i)<<std::endl;
  //   std::cout<<"erfc: " <<erfc(i)<<" E: "<<E(i)<<std::endl;
  // }

  // std::cout<<isSquare(1)<<sqrt(1)<<std::endl;

  return 0;
}




//y = 2000

//startTime = time.clock()
//for i in range(y):
//  x =  class_number(3, log(sqrt(3) + 2))
//print ((time.clock() - startTime)/y)
//print x

// cdef dict stored, units
// //started 1:14
// //f = open("fundamental_units_(first_half).txt")
// //f = open("medium_list")
// f = open("long_list")
// //f = open("short_list")
// //startTime = time.clock();
// answers = open("cython_list", "w")
// cdef list unit_list = str.split(f.read(), "\n")
// cdef list squares = [pow(x, 2) for x in range(1000)]
// cdef list components

// for s in unit_list:
//   components = str.split(s, " ", 1)
//   //units[<int>int(components[0])] = <unicode>(components[1])
//   units["a"] = "b"

// print "made list"
// print units

//cdef int startTime = time.clock()
//cdef int answer
//for d in units:
//  if d not in squares:
//      disc = d
//      for i in range(sqrt(disc)+1, 1, -1):
//        if disc % (i * i) == 0:
//          disc //= (i * i)

//      if disc == d:
//        answer = <int> class_number(d, log(units[d]))
//        units[d] = answer
//        answers.write(str(d) + " " + answer + "\n")
//        print d, answer
//  //print(str(start + step*(j+1)) + " finished")
//answers.close()
//f.close()
//print str(time.clock() - startTime) + " seconds"

//step = 100
//start = 300000
//startTime = time.clock();
//for j in range(1000):
//  solns = ""
//  for i in range(start + step*j, start + step*(j+1)):
//    if sqrt(i) not in ZZ:
//      //K = QuadraticField(i, 'x')
//      calculated = int(ty_sum(i, units, stored))
//      //if calculated != K.class_number():
//      //  print str(i) + " wrong "
//      //  solns += (str(i) + "\n")
//      //else:
//        //print str(i)
//      solns += (str(i) + " " + str(calculated)) + "\n"
//  answers.write(solns)
//  answers.flush()  
//  print(str(start + step*(j+1)) + " finished")
//answers.close()
//f.close()
//print str(time.clock() - startTime) + " seconds"