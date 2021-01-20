#include "FortranUtilsFC.h"
#include <chrono>
#include <thread>

void fu_sleep_microseconds(const int *duration)
{
   std::this_thread::sleep_for(std::chrono::microseconds(*duration));
}
