### Реализация подстановки предкомпилированных/динамически созданных типов в рантайме
#### типы передаются при компиляции + при запуске, например:
```bash
g++ -DDEBUG="true" -DTYPES="FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT" -DSIZES="S(1920,1080),S(36,84)" -std=c++2b fluid.cpp

./a.out --p-type="FLOAT" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FIXED(20, 10)" --n=36 --m=84 --field-file="field"
```

(течение жидкости написано не мной а преподавателем - моё только надстройки работы с типами Fixed, Fast_fixed, и размерами поля)

сравнение запусков с разным количеством потоков на разных по размеру полях:
![image](https://github.com/user-attachments/assets/95819030-2288-4783-b86f-a812c4225121)

базовый, не видоизменённый код на поле из примера 34 на 86 работает 1449.26 sec - смысла тестировать на больших полях я не вижу, тк это займёт ооочень много времени

Все тесты можно проверить использую команды в makefile:
```bash
make testBaseCode

make CreateBuildLib

make BuildRunCmakeRelease

make BuildRunCmakeReleaseSemiLarge

make BuildRunCmakeReleaseLarge
```