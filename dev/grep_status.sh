#!/bin/bash
cd /home/toki/toga2/TOGA2/src/python/modules
echo "== constants 330-380"; sed -n 320,380p cesar_wrapper_constants.py
echo "== conservation_summary 90-240"; sed -n 90,240p conservation_summary.py
echo "== input_producer N"; grep -n -B3 -A8 "\bN\b" input_producer.py | head -80
