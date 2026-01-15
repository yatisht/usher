#!/bin/sh
sed -i.bak 's/const target& this_ = \*this;/const class target\& this_ = *this;/' mutation_detailed.pb.cc
