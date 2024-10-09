#!/bin/bash
cd jbwa/jbwa-1.0.0/ && make && cd ../../ && sudo mvn compile assembly:single
cd ./target && sudo ln -s $(ls -t | head -n1) bisulfitehic-default.jar

if [ -z $BISHIC ]
then
	export BISHIC=${PWD}
	echo "export BISHIC=${BISHIC}" >> ~/.bash_profile
	source ~/.bash_profile
fi

cd .. && chmod 755 bisulfitehicMap