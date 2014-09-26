if [ "$1" = "clean" ]
then
	rm -rf build
else
	mkdir build
	cd build && cmake ..
	make -j 4
	cd ..
fi
if [ "$1" = "install" ]
then
	cp biphy /usr/local/bin/
fi
