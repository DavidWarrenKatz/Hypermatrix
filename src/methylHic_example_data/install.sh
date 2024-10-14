#!/bin/bash

# Step into the jbwa directory and run make, followed by Maven commands
cd jbwa/jbwa-1.0.0/ && make && cd ../../ && mvn compile assembly:single

# Set the BISHIC environment variable if not already set
if [ -z "$BISHIC" ]; then
    export BISHIC=${PWD}
    echo "export BISHIC=${BISHIC}" >> ~/.bash_profile
    source ~/.bash_profile
fi

# Ensure we're in the correct directory to create the symbolic link
cd ./target

# Check if the symbolic link already exists, and remove it if it does
if [ -L "bisulfitehic-default.jar" ]; then
    echo "Removing existing symbolic link 'bisulfitehic-default.jar'"
    rm bisulfitehic-default.jar
fi

# Identify the correct JAR file for linking (in case more files are in the directory)
jar_file=$(ls -t *.jar | head -n1)

# Ensure the JAR file was found
if [ -z "$jar_file" ]; then
    echo "Error: No JAR file found in the current directory."
    exit 1
fi

# Create the symbolic link
ln -s "$jar_file" bisulfitehic-default.jar

# Ensure the bisulfitehicMap script is executable
chmod 755 ../bisulfitehicMap

