#include "text_parser.hpp"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdio.h>

text_parser::text_parser(const std::string &path) {
    auto fp = fopen(path.c_str(), "rb");
    if (fp == NULL) {
        printf("File does not exist or cannot be opened.\n");
        exit(1);
    }

    fseek(fp, 0L, SEEK_END);
    size = ftell(fp); // text_parser size member
    rewind(fp);

    // Load all data into memory in one allocation
    data = new char[size];
    (void)fread((void *)data, 1, size, fp);
    fclose(fp);

    // Position iterators at correct locations
    curr = data;
    end = data + size;
    eol = curr;
    while (eol != end && *eol != '\n') {
        ++eol;
    }
}

text_parser::~text_parser() { delete[] data; }

void text_parser::next_line() {
    // Check end of file reached
    if (eol == end) {
        return;
    }
    // Otherwise, move current iterator one past end of line
    curr = eol + 1;

    // Move end of line iterator to current line eol
    eol = curr;
    while (eol != end && *eol != '\n') {
        ++eol;
    }
}

