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

# if $DEV, run a simpler build task
# else build the entire project
if [[ "$DEV" == "true" ]]; then
	cd "test/build"
		if build; then
			echo
			./geometry_test
		fi
	cd ../../
else
	cd "test/build"
		build
	cd ../../
fi