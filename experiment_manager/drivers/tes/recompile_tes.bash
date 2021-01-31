#!/bin/bash

(cd ./../../../src && make clean)
(cd ./../../../src && make)
cp ./../../../src/output ./
