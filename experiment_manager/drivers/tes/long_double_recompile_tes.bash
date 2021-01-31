#!/bin/bash


(cd ./../../../src && make clean)
(cd ./../../../src && make long)
cp ./../../../src/output ./output_longdouble
