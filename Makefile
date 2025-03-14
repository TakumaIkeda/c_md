# コンパイラ
CC = gcc

# コンパイラオプション
CFLAGS = -Wall -I./src

# 実行ファイル名
TARGET = md

# ソースファイルとオブジェクトファイルのリスト
SRCS = main.c $(wildcard src/*.c)
OBJS = $(SRCS:.c=.o)

# 実行ファイルの生成ルール
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm -O3

# オブジェクトファイルの生成ルール
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# クリーンルール
clean:
	rm -f $(TARGET) $(OBJS)