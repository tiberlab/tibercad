
TARGET1 = drawlines
TARGET2 = gen2xyz
TARGET = $(TARGET1) $(TARGET2)

all: $(TARGET) 

$(TARGET1): 
	gcc -o $@ drawlines.c -lm

$(TARGET2): 
	gcc -o $@ gen2xyz.c -lm

clean:
	rm -f $(TARGET)

