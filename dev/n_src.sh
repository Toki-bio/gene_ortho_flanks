#!/bin/bash
cd /home/toki/toga2/TOGA2/src/python/modules
sed -n 800,870p input_producer.py
grep -n "ILLEGAL\|illegal" input_producer.py | head
