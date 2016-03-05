#include <cmath>
#include "gaussianCDF.h"

double gaussianCDF(double x) {
  return 0.5 * std::erfc(-x * M_SQRT1_2);
}