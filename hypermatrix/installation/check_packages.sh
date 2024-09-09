#If you are having problems installing the enviroment on your computers, 
#try running this file.

#!/bin/bash

# Check if hic-straw and fanc are installed
echo "Checking if hic-straw and fanc are installed..."

if ! pip list | grep -q hic-straw; then
  echo "hic-straw not installed. Installing with --use-pep517..."
  pip install --use-pep517 hic-straw || { echo "Failed to install hic-straw"; exit 1; }
fi

if ! pip list | grep -q fanc; then
  echo "fanc not installed. Installing with --use-pep517..."
  pip install --use-pep517 fanc || { echo "Failed to install fanc"; exit 1; }
fi


