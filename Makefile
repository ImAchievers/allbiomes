CC      = gcc
CFLAGS  = -O3 -march=native -fwrapv -flto -Wall -Ilib
LDFLAGS = -lm -lpthread

CORE_SRCS = lib/generator.c lib/layers.c lib/biomenoise.c lib/biomes.c lib/noise.c lib/terrainnoise.c
CORE_OBJS = $(CORE_SRCS:.c=.o)

all: allbiomes

allbiomes: allbiomes.c $(CORE_OBJS)
	$(CC) $(CFLAGS) allbiomes.c $(CORE_OBJS) -o allbiomes $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

run: allbiomes
	./allbiomes

clean:
	rm -f allbiomes allbiomes.exe $(CORE_OBJS)

.PHONY: all run clean
