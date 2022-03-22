trap exit SIGINT

while true;
	do find src/ | entr -d sh ./bin/build.sh -d; 
done