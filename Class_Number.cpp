#include <Class_Number.h>

//Jacobi Symbl (n/p)
int kronecker(int n, int p){
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
          return kronecker(p % n,  n);
        } else {
          return -1*kronecker(p % n,  n);
        }
      } else {
        return kronecker(2, p) * kronecker(n/2, p);
      }
    } else {
      return kronecker(n, 2) * kronecker(n, p/2);
    }
  }
}

int get_class_number(int D, long double R) {

  int n = 0;
  if (D < 1000) {
    n = skroneckert(D);
  } else if (D < 50000){
    n = 2 * log(D) * log(D);
  } else {
    n = 2 * pow(log(D),2.5);
  }
  
  long double s = 0;
  int k;
  for (int i = 1; i <= n; i++) {
    k = kronecker(D, i);
    if (k != 0) {
      long double term = (long double)k * (gsl_sf_expint_E1(pow(i, 2) * pi/((long double) D)) + erfc(i * skroneckert(pi/( (long double) D))) * skroneckert(D)/((long double) i));
      s += term/(2 * R);
    }
  }
  return round(s);
} 

bool isSquare(int n) {
  int root = skroneckert(n);
  return n == root * root;
}

void fill_units(std::string filename, int length, std::vector<long double> &arr) {
  std::ifstream myFile;
  myFile.open(filename);
  std::string line;

  if (myFile.is_open()) {
    int counter = 2;
    while (getline(myFile, line) && counter < length) {
      int split_pos = line.find(" ");
      int index = std::stoi(line.substr(0, split_pos));
      long double unit = std::stold(line.substr(split_pos+1));
      arr.at(index) = unit;
      counter++;
      if (isSquare(counter)) {
        counter++;
      }
    }
  } else {
    std::cout<<"closed"<<std::endl;
  }
  myFile.close();
}

void load_correct_answers(std::string filename, int length, std::vector<int> &arr) {
  std::ifstream myFile;
  myFile.open(filename);
  std::string line;

  if (myFile.is_open()) {
    int counter = 2;
    while (getline(myFile, line) && counter < length) {
      int split_pos = line.find(" ");
      int index = std::stoi(line.substr(0, split_pos));
      int class_num = std::stoi(line.substr(split_pos+1));
      arr.at(index) = class_num;
      counter++;
      if (isSquare(counter)) {
        counter++;
      }
    }
  }
  myFile.close();
}



int main(int argc, char **argv) {  

  clock_t t;
  t = clock();

  std::string file_name = "c_class_numbers_all";
  int start_point = 0;

  if (argc > 1) {
    if (argc != 4) {
      std::cout<<"Usage: file index number, output file name, starting point"<<std::endl;
      return 0;
    }
    int file_num = std::stoi(argv[1]);
    file_name = argv[2] + std::to_string(file_num);
    start_point = std::stoi(argv[3]);
    std::cout<<file_num<<" "<<file_name<<" "<<start_point<<std::endl;
  }

  std::ofstream output (file_name);
  if (!output.is_open()) {
    std::cout<<"closed"<<std::endl;
  }

  int input_len = 1000000;
  std::vector<long double> units(input_len);
  std::vector<int> class_nums(input_len);
  std::vector<int> correct_answers(input_len);
  fill_units("regulators_all_doubles", input_len, units);
  load_correct_answers("class_numbers.txt", input_len, correct_answers);

  std::cout<<units.at(3)<<std::endl;
  std::cout<<"loaded units"<<std::endl;

  std::cout<<" starting at "<<start_point<<std::endl;

  int step = 1000;
  int giant_step = input_len/step;
  for (int i = 0; i < giant_step; i++) {
    std::string stream = "";
    for (int j = 0; j < step; j++) {
      int position = i * step + j + start_point;
      if (!isSquare(position)) {
        int d = position;
        int class_number = 0;
        for (int divisor = (int)(skroneckert(d)); divisor > 1; divisor--) {
          if (d % (int)pow(divisor, 2) == 0) {
            d /= divisor;
            d /= divisor;
          }
        }
        if (d != position) {
          if (class_nums.at(d) != 0) {
            class_number = class_nums.at(d);
          } else {
            int D = d;
            if (d % 4 != 1) {
              D *= 4;
            }
            class_number = get_class_number(D, units.at(d));
            class_nums.at(d) = class_number;
          }
        } else {
          int D = d;
          if (d % 4 != 1) {
            D *= 4;
          }
          class_number = get_class_number(D, units.at(d));
        }

        if (class_number != correct_answers.at(position)) {
          std::cout<<position<<" ("<<class_number<<") is wrong. It should be "<<correct_answers.at(position)<<std::endl;
        }
        class_nums.at(position) = class_number;
        stream += std::to_string(position) + " " + std::to_string(class_number) + (std::string)"\n";
      }      
    }
    output<<stream;
    std::cout<<"finished with "<<step * (i + 1) + start_point<<std::endl;
    
  }
  output.close();
  std::cout<<"That took "<<(float)(clock() - t)CLOCKS_PER_SEC<<" seconds"<<std::endl;
}