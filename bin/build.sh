build() {
	cmake -S test -B test/build
	cmake --build test/build
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
if [ "$DEV" == "true" ] && build; then
	echo
	./test/build/geometry_test
else
	build
fi