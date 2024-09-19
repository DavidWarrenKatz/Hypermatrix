#!/bin/bash

# Variables for the paths
BISULFITEHIC_DIR="/your/path/to/bisulfitehic"
LIBRARY_PATH="$BISULFITEHIC_DIR/jbwa/src/main/native"
JAR_PATH="$BISULFITEHIC_DIR/target/bisulfitehic-0.38-jar-with-dependencies.jar:$BISULFITEHIC_DIR/lib/jbwa.jar"

# Run the help command for bisulfite Hi-C mapping
echo "Running bisulfite Hi-C help command..."

java -Xmx15G -Djava.library.path="$LIBRARY_PATH" \
-cp "$JAR_PATH" \
main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem --help

