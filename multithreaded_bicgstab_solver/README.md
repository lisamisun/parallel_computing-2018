# Многопоточная реализация солвера BiCGSTAB для СЛАУ с разреженной матрицей, заданной в формате CSR

Текст задания находится в файле *task1.pdf*, отчет по его выполнению находится в файле *report.pdf*.

Программа написана на языке ```C++```.

   Программа собирается командой: **make main**  
   Запустить программу: **./main (+ работают параметры командной строки, описанные в задании)**  
   Собрать тесты корректности базовых операций: **make test**  
   Запустить тесты корректности базовых операций: **./test**  
   Удалить main и test: **make clean**  
   
   Чтобы файл test скомпилировался, надо скачать ```googletest```:  
   **git clone https://github.com/google/googletest.git**  
   **cd googletest**  
   **cmake .**  
   **make**  

___

# Multithreaded implementation of the BiCGSTAB solver for a linear system with a sparse matrix specified in CSR format

The text of the task is in the *task1.pdf* file and its implementation report is in the *report.pdf* file.
