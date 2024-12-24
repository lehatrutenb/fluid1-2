Runtime set type/sz to use precompile/dinamic created objects

(течение жидкости написано не мной - моё только надстройки работы с типами Fixed, Fast_fixed, и размерами)

сравнение запусков с разным количеством потоков на разных по размеру полях:
![image](https://github.com/user-attachments/assets/95819030-2288-4783-b86f-a812c4225121)

базовый, не видоизменённый код на поле из примера 34 на 86 работает 1449.26 sec - смысла тестировать на больших полях я не вижу, тк это займёт ооочень много времени

Все тесты можно проверить использую команды в makefile:
make testBaseCode
make CreateBuildLib
make BuildRunCmakeRelease
make BuildRunCmakeReleaseSemiLarge
make BuildRunCmakeReleaseLarge
