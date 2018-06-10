CC = gcc

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

nn3d: jsmn.o pq.o kdtree.o nn3d.o
	$(CC) jsmn.o pq.o kdtree.o nn3d.o -o nn3d

.PHONY: clean

clean:
	rm *.o nn3d
