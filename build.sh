if [ "$1" = "clean" ]; then
	rm -rf build
else
if [ ! -d build ]; then
	mkdir build
fi	
	cd build && cmake ..
	make -j 4
	cd ..
fi
if [ "$1" = "install" ]
then
	cp biphy /usr/local/bin/
fi
