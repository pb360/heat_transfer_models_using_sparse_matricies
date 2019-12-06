#include <iostream>
#include <string>

#include "heat.hpp"
#include "sparse.hpp"




int main(int argc, char *argv[])
{
  /* Get command line arguments */
  if (argc != 3)
  {
    std::cout << "Usage:" << std::endl;
    std::cout << "  " << argv[0] << " <input file> <soln prefix>" << std::endl;
    return 0;
  }
  
  std::string inputfile     = argv[1];
  std::string soln_prefix   = argv[2];

  SparseMatrix CSR(inputfile, soln_prefix);
  

  std::cout << "main ran fully \n \n"; 

  return 0;
}
       