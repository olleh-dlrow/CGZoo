/*
clang++ getch.cpp -lpthread -o getch
*/

#include<stdio.h>
#include<unistd.h>
#include<pthread.h>
#include<fcntl.h>
#include <termios.h>

int _getch(void) {
    int cr;
    struct termios nts, ots;

    if (tcgetattr(0, &ots) < 0) // 得到当前终端(0表示标准输入)的设置
        return EOF;

    nts = ots;
    cfmakeraw(&nts); // 设置终端为Raw原始模式，该模式下所有的输入数据以字节为单位被处理
    if (tcsetattr(0, TCSANOW, &nts) < 0) // 设置上更改之后的设置
        return EOF;

    cr = getchar();
    if (tcsetattr(0, TCSANOW, &ots) < 0) // 设置还原成老的模式
        return EOF;

    return cr;
}

void* print(void* args) {
    int* pcol = (int*)args;
    printf("color: %d\n", *pcol);
    char ch;
    do {
        ch = _getch();
        switch (ch)
        {
        case 'a':
            *pcol -= 5;
            break;
        case 'd':
            *pcol += 5;
            break;
        default:
            break;
        }
        printf("color: %d\n", *pcol);
    }while(ch != 'q');
    return 0;
}

int main() {
    pthread_t tid;
    pthread_attr_t attr;
    int color = 0;

    pthread_attr_init(&attr);
    pthread_create(&tid, &attr, print, &color);
    pthread_join(tid, NULL);
    return 0;
}