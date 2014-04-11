#pragma once
#ifdef  _WIN32
#define getcwd(x,y) GetCurrentDirectory((y),(x))
#else
#include <unistd.h>
#include <arpa/inet.h>
#endif 


