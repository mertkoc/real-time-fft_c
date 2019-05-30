

output:
	gcc -shared -o overlapsave.so -fPIC overlapsave.c

clean: 
	rm *.so
