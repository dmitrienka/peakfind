# README

## О пргограмме

peakfind.r - R-скрипт, облегчающий подготовку порошковых CIF для ICDD. Использует готовое уточнение в TOPAS4-2.

## Использование

Обработайте дифрактограмму в TOPAS4-2. Целевая фаза должна быть описана как _hkl Phase_. Убедитесь, что режимы _Show calculated curve_, _Show difference curve_ и _Show background curve_ активны (после включения режимов необходимо выполнить _Refine_). Экспортируйте данные в текстовый файл (правой кнопкой мыши по паттерну и _Save if displayed..._), например **data.txt**. Сохраните проект, например **icdd.pro**.

Вызывайте скрипт следующим образом:
 
 ````
 > peakfind.r data.txt icdd.pro
 ````
 Будут созданы два файла: диагностическая картинка **check.eps** и каркас CIF **icdd.cif**.
 
## Принцип действия
 
 * Скрипт находит максимумы в расчетных данных
 * В окресности максимумов расчетных данных находятся максимумы экспериментальных данных
 * 
