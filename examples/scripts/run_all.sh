#!/bin/bash

for i in ex_*.sh; do
  sbatch ${i}
done
