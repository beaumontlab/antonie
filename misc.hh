#pragma once
#include <stdio.h>
#include <string>
#include <stdint.h>

void chomp(char* line);
char* sfgets(char* p, int num, FILE* fp);
void reverseNucleotides(std::string* nucleotides);
uint64_t filesize(const char* name);
bool stringfgets(FILE* fp, std::string* line);
