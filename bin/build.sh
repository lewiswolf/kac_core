build() {
	cmake -S ../ -B .
	make
}

# loop over args and check for -d
for arg in "$@"
do
	case "$arg" in
		-d|--dev)
			DEV="true"
			;;
		*)
			;;
	esac
done

# build tests and run if in dev mode 
cd "test/build"
	if [ "$DEV" == "true" ] && build; then
		echo
		./geometry_test
	else
		build
	fi
cd ../../