#Makefile
#Edit for your compiler


all: 
	mkdir -p output
	cd src && make && cd ..

clean:
	cd src && make clean && cd ..
cleaner:
	cd src && make cleaner && cd ..
